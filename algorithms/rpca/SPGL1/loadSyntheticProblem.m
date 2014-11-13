function [Y,L,S,params] = loadSyntheticProblem(varargin)
% [Y,L,S,params] = loadSyntheticProblem(problemName, referenceSolutionDirectory)
%
% Returns parameters and solution for the following problems
%   (parameters are tuned so that each problem has the
%    the same true solution (L,S) )
% Or, call this with no input arguments to set the path...
%
% "Lag": min   lambdaL*||L||_* + lambdaS*||S||_1 + .5||L+S-Y||_F^2
% 
% "Sum": min   ||L||_* + lambdaSum*||S||_1 
%        subject to ||L+S-Y||_F <= epsilon
%
% "Max": min   max( ||L||_* , lambdaMax*||S||_1  )
%        subject to ||L+S-Y||_F <= epsilon
%
% "Dual-sum"  min  .5||L+S-Y||_F^2
%         subject to ||L||_* + lambdaSum*||S||_1  <= tauSum
%
% "Max-sum"  min  .5||L+S-Y||_F^2
%         subject to max( ||L||_*, lambdaMax*||S||_1)  <= tauMax
%
% so params has the following fields:
%   .lambdaL
%   .lambdaS
%   .lambdaSum
%   .lambdaMax
%   .tauSum
%   .tauMax
%   .epsilon
%
%   .objective  This is a function of (L,S) and is the objective
%               of the "Lagrangian" formulation
%
% Stephen Becker, March 6 2014



% addpath for the proxes and solvers
Mfilename = mfilename('fullpath');
baseDirectory = fileparts(Mfilename);
addpath(genpath(baseDirectory));

if nargin==0
    % we just setup the path, nothing else...
    return;
end

if nargin >= 2
    baseDirectory = varargin{2};
else
    baseDirectory = fullfile(baseDirectory,'referenceSolutions');
end
norm_nuke = @(L) sum(svd(L,'econ'));
vec = @(X) X(:);
problemString = varargin{1};

%%
% -- Common problem strings --
% problemString = 'small_randn_v2';
% problemString = 'medium_randn_v1';
% problemString = 'medium_exponentialNoise_v1';
% problemString = 'large_exponentialNoise_square_v1';
% problemString = 'large_exponentialNoise_rectangular_v1';
% problemString = 'huge_exponentialNoise_rectangular_v1'; % not yet done!
% problemString = 'nsa500_SNR80_problem1_run1';

fprintf('%s\n** Problem %s **\n\n',char('*'*ones(60,1)), problemString );
fileName      = fullfile(baseDirectory,[problemString,'.mat']);
if exist(fileName,'file')
    fprintf('Loading reference solution from file\n');
    load(fileName); % should contain L, S and Y at least, and preferably lambdaL and lambdaS
else
    [L,S,Y] = deal( [] );
end

% For problems from the NSA software package:
if length(problemString)>3 && strcmpi( problemString(1:3), 'nsa' )
    % This sscanf line is broken in versions of Matlab after R2010a
    tmp = sscanf(problemString,'nsa%d_SNR%d_problem%d_run%d');
    m       = tmp(1); % either 500, 1000 or 1500
    SNR     = tmp(2); % either 80 or 45
    problem_type = tmp(3);
    run_counter  = tmp(4);
    % e.g.,  'nsa500_SNR80_problem1_run1'
    n  = m;
%     % Variable parameters
%     m = 500; n = m;
%     problem_type = 1; % between 1 and 4. Controls rank/sparsity
%     SNR          = 80;
%     run_counter  = 1; % between 1 and 10. Only affects seed
            
            
%     % Possible sparsity amounts:
%     p_array = ceil( m*n*[ 0.05, 0.1, 0.05, 0.1 ] );
%     r_array = ceil(min(m,n)*[  0.05, 0.05, 0.1, 0.1 ] );
    
    % Changing the conventions, March 14. They had 4 cases, for p = {.05,.1}
    %   and r = {.05, .1 } (these are sparsity and rank percentages)
    % I don't want to do that many tests since we will display graphically
    %   and not in tables. I also want more values of p and r. So I will
    %   make p and r be linked.
    p_array = ceil( m*n*[ 0.05, 0.1, 0.2, 0.3 ] );
    r_array = ceil(min(m,n)*[  0.05, 0.1, 0.2, 0.3 ] );
    
    p = p_array(problem_type);
    r = r_array(problem_type);
    pwr = r+p/(m*n)*100^2/3;
    stdev = sqrt(pwr/10^(SNR/10));
    delta = sqrt(min(m,n)+sqrt(8*min(m,n)))*stdev;
    lambdaSum = 1/sqrt( max(m,n) ); % "xi" in NSA
else
    switch problemString
        case 'small_randn_v1'
            my_rng(101);
            m   = 40; n = 42; lambdaL = .5; lambdaS = 1e-1;
            if isempty(Y)
                Y   = randn(m,n);
                printEvery = 500; tol = 1e-12;  maxIts = 3e3;
            end
            
        case 'small_randn_v2'
            my_rng(102);
            m   = 40; n = 42; lambdaL = .5; lambdaS = 1e-1;
            if isempty(Y)
                Y   = randn(m,n);
                printEvery = 500; tol = 1e-12; maxIts = 3e3;
            end
        case 'medium_randn_v1'
            my_rng(101);
            m   = 400; n = 500; lambdaL = .25; lambdaS = 1e-2;
            if isempty(Y)
                Y   = randn(m,n);
                printEvery = 1; %tol = 1e-4; maxIts = 100;
                tol = 1e-10; maxIts = 1e3;
            end
        case 'medium_exponentialNoise_v1'
            % Converges much faster than 'medium_randn_v1'
            my_rng(101);
            m   = 400; n = 500; lambdaL = .25; lambdaS = 1e-2;
            rk  = round(.5*min(m,n) );
            Q1  = haar_rankR(m,rk,false);
            Q2  = haar_rankR(n,rk,false);
            Y   = Q1*diag(.1+rand(rk,1))*Q2';
            mdn = median(abs(Y(:)));
            Y   = Y + exprnd( .1*mdn, m, n );
            
            printEvery = 5; tol = 1e-10; maxIts = 1000;
            
        case 'medium_exponentialNoise_v2'
            % Converges much faster than 'medium_randn_v1'
            % our quasi-Newton scheme does well here
            my_rng(101);
            m   = 400; n = 500; lambdaL = .25; lambdaS = 1e-2;
            rk  = round(.05*min(m,n) );
            Q1  = haar_rankR(m,rk,false);
            Q2  = haar_rankR(n,rk,false);
            Y   = Q1*diag(.1+rand(rk,1))*Q2';
            mdn = median(abs(Y(:)));
            Y   = Y + exprnd( .1*mdn, m, n );
            
            printEvery = 5; tol = 1e-10; maxIts = 1000;
            
        case 'medium_exponentialNoise_v3'
            % Converges much faster than 'medium_randn_v1'
            my_rng(101);
            m   = 400; n = 500; lambdaL = .25; lambdaS = 1e-2;
            rk  = round(.1*min(m,n) );
            Q1  = haar_rankR(m,rk,false);
            Q2  = haar_rankR(n,rk,false);
            Y   = Q1*diag(.1+rand(rk,1))*Q2';
            mdn = median(abs(Y(:)));
            Y   = Y + exprnd( .1*mdn, m, n );
            
            printEvery = 5; tol = 1e-10; maxIts = 1000;
            
        case 'large_exponentialNoise_square_v1'
            my_rng(201);
            m   = 1e3; n = m; lambdaL = .25; lambdaS = 1e-2;
            rk  = round(.5*min(m,n) );
            Q1  = haar_rankR(m,rk,false);
            Q2  = haar_rankR(n,rk,false);
            Y   = Q1*diag(.1+rand(rk,1))*Q2';
            mdn = median(abs(Y(:)));
            Y   = Y + exprnd( .1*mdn, m, n );
            
            printEvery = 5; tol = 1e-8; maxIts = 250;
            
        case 'large_exponentialNoise_rectangular_v1'
            my_rng(201);
            m   = 5e3; n = 200; lambdaL = .25; lambdaS = 5e-3; % OK
            rk  = round(.5*min(m,n) );
            Q1  = haar_rankR(m,rk,false);
            Q2  = haar_rankR(n,rk,false);
            Y   = Q1*diag(.1+rand(rk,1))*Q2';
            mdn = median(abs(Y(:)));
            Y   = Y + exprnd( .1*mdn, m, n );
            
            printEvery = 5; tol = 1e-8; maxIts = 100;
            
        case 'huge_exponentialNoise_rectangular_v1';
            disp('Warning: lambda values not yet tuned');
            my_rng(301);
            m   = 5e4; n = 200; lambdaL = .25; lambdaS = 5e-2;
            rk  = round(.5*min(m,n) );
            Q1  = haar_rankR(m,rk,false);
            Q2  = haar_rankR(n,rk,false);
            Y   = Q1*diag(.1+rand(rk,1))*Q2';
            mdn = median(abs(Y(:)));
            Y   = Y + exprnd( .1*mdn, m, n );
            
            printEvery = 5; tol = 1e-10; maxIts = 100;
            
        
            
        otherwise
            error('bad problem name');
    end
end


% Generate new reference soln if not precomputed
if isempty(L)
    fprintf('Generating reference solution...\n');
    if length(problemString)>3 && strcmpi( problemString(1:3), 'nsa' )
        % We solve Sum problem, which gives us a lot of parameters but not
        % the parameters for the Lagrangian problem unfortunately.
        
        % I am changing the seeds format
        seed = 601;
        my_rng( seed + run_counter );
        [Y, optimal_L, optimal_S, deltabar]=create_data_noisy_L2(m,n,p,r,stdev);
        
        tol = 0.0001;
        denoise_flag=0;

        %         [L,S,out] = nsa(Y,stdev,tol,denoise_flag,optimal_L,optimal_S);
        %         [L,S,out] = nsa(Y,-delta,tol,denoise_flag,optimal_L,optimal_S);
       
        [L,S,out] = nsa(Y,-delta,tol,denoise_flag,[],[],lambdaSum);
        fprintf('Discrepancy with intitial signal L is %.2e\n', ...
            norm(L-optimal_L,'fro')/norm(L,'fro') );
        fprintf('Discrepancy with intitial signal S is %.2e\n', ...
            norm(S-optimal_S,'fro')/norm(S,'fro') );
        
        lambdaL     = [];
        lambdaS     = [];
        maxIts      = [];
    else
        tic
        opts    = struct('printEvery',printEvery,'FISTA',true,'tol',tol,'maxIts',maxIts);
        %[L,S,errHist] = solver_RPCA_Lagrangian_simple(Y,lambdaL,lambdaS,[],opts);
        opts.quasiNewton  = true; opts.FISTA=false; opts.BB = false;
        otps.SVDstyle   = 1; % for highest accuracy
        [L,S,errHist] = solver_RPCA_Lagrangian(Y,lambdaL,lambdaS,[],opts);
        toc
    end
    
    save(fileName,'L','S','Y','lambdaL','lambdaS','tol','maxIts');
end

sigma   = svd(L,'econ');
rk      = sum( sigma > 1e-8 );
normL   = sum(sigma);
normS   = norm(S(:),1);

params = [];
params.objective    = @(L,S) lambdaL*norm_nuke(L)+lambdaS*norm(S(:),1) + ...
    .5*norm(vec(L+S-Y))^2;
params.epsilon     = norm( L + S - Y, 'fro' );
params.lambdaL     = lambdaL;
params.lambdaS     = lambdaS;
if exist('lambdaSum','var')
    params.lambdaSum   = lambdaSum;
else
    params.lambdaSum   = lambdaS/lambdaL;
end
params.lambdaMax   = normL/normS;
params.tauSum      = normL + params.lambdaSum*normS;
params.tauMax      = max(normL,params.lambdaMax*normS );

fprintf('\tSize %d x %d\n\tL has nuclear norm %.3f and rank %d of %d possible\n', ...
    m,n,normL, rk, length(sigma) );
fprintf('\tS has l1 norm %.3f and %.2f%% of its elements are nonzero\n', ...
    normS, 100*sum( abs(S(:)) > 1e-8)/numel(S) );
fprintf('\t||L+S-Y||/||Y|| is %.2e\n', params.epsilon/norm(Y,'fro' ) );



end % end of main function


function Q = haar_rankR(n, r, cmplx)
% Q = HAAR_RANKR( N , R)
%   returns a N x N unitary matrix Q
%   from the Haar measure on the
%   Circular Unitary Ensemble (CUE)
%
% Q = HAAR_RANKR( N, R, COMPLEX )
%   if COMPLEX is false, then returns a real matrix from the
%   Circular Orthogonal Ensemble (COE). By default, COMPLEX = true
% 
% For more info, see http://www.ams.org/notices/200705/fea-mezzadri-web.pdf
% "How to Generate Random Matrices fromthe Classical Compact Groups"
% by Francesco Mezzadri (also at http://arxiv.org/abs/math-ph/0609050 )
%
% Other references:
%   "How to generate a random unitary matrix" by Maris Ozols
%
% Stephen Becker, 2/24/11. updated 11/1/12 for rank r version.
%
% To test that the eigenvalues are evenly distributed on the unit
% circle, do this:   a = angle(eig(haar(1000))); hist(a,100);
%
% This calls "randn", so it will change the stream of the prng

if nargin < 3, cmplx = true; end
if nargin < 2 || isempty(r) , r=n; end

if cmplx
    z = (randn(n,r) + 1i*randn(n,r))/sqrt(2.0);
else
    z = randn(n,r);
end
[Q,R] = qr(z,0); % time consuming
d = diag(R);
ph = d./abs(d);
% Q = multiply(q,ph,q) % in python. this is q <-- multiply(q,ph)
%   where we multiply each column of q by an element in ph

%Q = Q.*repmat( ph', n, 1 );
%  use bsxfun for this to make it fast...
Q = bsxfun( @times, Q, ph' );

end % end of haar_rankR


function varargout = my_rng(varargin)
%MY_RNG Control the random number generator used by RAND, RANDI, and RANDN (SRB version)
%   MY_RNG(SD) seeds the random number generator using the non-negative
%   integer SD so that RAND, RANDI, and RANDN produce a predictable
%   sequence of numbers.
% MODIFIED BY STEPHEN BECKER
%   See also RAND, RANDI, RANDN, RandStream, NOW.


%   See <a href="http://www.mathworks.com/access/helpdesk/help/techdoc/math/brn4ixh.html#brvku_2">Choosing a Random Number Generator</a> for details on these generators.

if exist('rng','file')
    [varargout{1:nargout}] = rng(varargin{:});
    return;
end
% Otherwise, imitate the functionality...

% -- SRB adding this --
error(nargchk(1,1,nargin));
error(nargoutchk(0,0,nargout));
arg1 = varargin{1};
% For R2008a, this doesn't work... (not sure what earliest version is)
if verLessThan('matlab','7.7')
    randn('state',arg1);
    rand('state',arg1);
else
    if verLessThan('matlab','8.2')
        RandStream.setDefaultStream(RandStream('mt19937ar', 'seed', arg1) );
    else
        RandStream.setGlobalStream(RandStream('mt19937ar', 'seed', arg1) );
    end
end
end % end of my_rng
