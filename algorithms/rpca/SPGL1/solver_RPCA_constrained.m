function [L,S,errHist] = solver_RPCA_constrained(AY,lambda_S, tau, A_cell, opts)
% [L,S,errHist] = solver_RPCA_constrained(Y,lambda_S, tau, A_cell, opts)
% Solves the problem
%   minimize_{L,S} .5|| L + S - Y ||_F^2 
%   subject to 
%   if opts.sum = true
%       (1)  ||L||_* + lambda_S ||S||_1 <= tau
%   if opts.max = true
%       (2) max(  ||L||_* , lambda_S ||S||_1 ) <= tau
%
%   if opts.max and opts.sum are false and tau is a negative number, then
%   we solve the problem:
%       minimize_{L,S} .5|| L + S - Y ||_F^2  + abs(tau)*( ||L||_* + lambda_S ||S||_1 )
%   (but see solver_RPCA_Lagrangian.m for a simpler interface)
%
%   or if A_cell is provided, where A_cell = {A, At}
%   (A is a function handle, At is a function handle to the transpose of A)
%   then
%
%   minimize_{L,S} .5|| A(L + S) - Y ||_F^2 
%       subject to ...
%   (here, Y usually represents A(Y); if Y is not the same size
%    as A(L), then we will automatically set Y <-- A(Y) )
%
%   errHist(1,:) is a record of the residual
%   errHist(2,:) is a record of the full objective (that is, .5*resid^2 )
%   errHist(3,:) is the output of opts.errFcn if provided
%
% opts is a structure with options:
%   opts.sum, opts.max  (as described above)
%   opts.L0         initial guess for L (default is 0)
%   opts.S0         initial guess for S (default is 0)
%   opts.size       [n1,n2] where L and S are n1 x n2 matrices. The size is automatically
%       determined in most cases, but when providing a linear operator
%       it may be necessary to provide an explicit size.
%   opts.tol        sets stopping tolerance
%   opts.maxIts     sets maximum number of iterations
%   opts.printEvery will print information this many iterations
%   opts.displayTime will print out timing information (default is true for large problems)
%   opts.errFcn     a function of (L,S) that records information
%   opts.trueObj    if provided, this will be subtracted from errHist(2,:)
%   opts.Lip        Lipschitz constant, i.e., 2*spectralNorm(A)^2
%                       by default, assume 2 (e.g., good if A = P_Omega)
%   opts.FISTA      whether to use FISTA or not. By default, true
%     opts.restart  how often to restart FISTA; set to -Inf to make it automatic
%   opts.BB         whether to use the Barzilai-Borwein spectral steplength
%     opts.BB_type  which BB stepsize to take. Default is 1, the larger step
%     opts.BB_split whether to calculate stepslengths for S and L independently.
%       Default is false, which is recommended.
%   opts.quasiNewton  uses quasi-Newton-like Gauss-Seidel scheme.
%                     Only available in "max" mode
%     opts.quasiNewton_stepsize     stepsize length. Default is .8*(2/Lip)
%     opts.quasinewton_SLS          whether to take S-L-S sequence (default is true)
%                                   otherwise, takes a L-S Gauss-Seidel sequence
%   opts.SVDstyle   controls what type of SVD is performed.
%       1 = full svd using matlab's "svd". Best for small problems
%       2 = partial svd using matlab's "svds". Not recommended.
%       3 = partial svd using PROPACK, if installed. Better than option 2, worse than 4
%       4 = partial svd using randomized linear algebra, following
%           the Halko/Tropp/Martinnson "Structure in Randomness" paper
%       in option 4, there are additional options:
%       opts.SVDwarmstart   whether to "warm-start" the algorithm
%       opts.SVDnPower  number of power iterations (default is 2 unless warm start)
%       opts.SVDoffset  oversampling, e.g., "rho" in Tropp's paper. Default is 5
%
% Stephen Becker, March 6 2014. Edited March 14 2-14. stephen.beckr@gmail.com
% See also solver_RPCA_Lagrangian.m, solver_RPCA_SPGL1.m


% todo: allow S >= 0 constraints, since this is easy

%error(nargchk(3,5,nargin,'struct'));
if nargin < 5, opts = []; end
% == PROCESS OPTIONS ==
function out = setOpts( field, default )
    if ~isfield( opts, field )
        opts.(field)    = default;
    end
    out = opts.(field);
    opts    = rmfield( opts, field ); % so we can do a check later
end

if nargin < 4 || isempty(A_cell)
    A   = @(X) X(:);
    [n1,n2] = size(AY);
    At  = @(x) reshape(x,n1,n2);
    
    if ~iscell(tau)
        AY  = A(AY);
    end
    % The factor of 2 is since we are in both L and S 
else
    A   = A_cell{1};
    At  = A_cell{2};
    % Y could be either Y or A(Y)
    if size(AY,2) > 1
        % AY is a vector, so it is probably Y and not AY
        disp('Changing Y to A(Y)');
        AY = A(AY);
        [n1,n2] = size(AY); % April 24 '14
    else
        % April 24, '14: we need to know the (n1,n2)
        sz      = setOpts('size',[] );
        if isempty(sz)
            error('Cannot determine the size of the variables; please specify opts.size=[n1,n2]');
        end
        n1 = sz(1);
        n2 = sz(2);
    end
    %[n1,n2] = size(AY); % comment out April 24 '14
end

% Some problem sizes. Feel free to tweak. Mainly affect the defaults
SMALL   = ( n1*n2 <= 50^2 );
MEDIUM  = ( n1*n2 <= 200^2 ) && ~SMALL;
LARGE   = ( n1*n2 <= 1000^2 ) && ~SMALL && ~MEDIUM;
HUGE    = ( n1*n2 > 1000^2 );


% -- PROCESS OPTIONS -- (some defaults depend on problem size )
tol     = setOpts('tol',1e-6*(SMALL | MEDIUM) + 1e-4*LARGE + 1e-3*HUGE );
maxIts  = setOpts('maxIts', 1e3*(SMALL | MEDIUM ) + 400*LARGE + 200*HUGE );
printEvery  = setOpts('printEvery',100*SMALL + 50*MEDIUM + 5*LARGE + 1*HUGE);
errFcn      = setOpts('errFcn', [] );
Lip         = setOpts('Lip', 2 );
restart     = setOpts('restart',-Inf);
trueObj     = setOpts('trueObj',0);
sumProject  = setOpts('sum', false );
maxProject  = setOpts('max', false );
QUASINEWTON = setOpts('quasiNewton', maxProject );
FISTA       = setOpts('FISTA',~QUASINEWTON);
BB          = setOpts('BB',~QUASINEWTON);
BB_split    = setOpts('BB_split',false);
BB_type     = setOpts('BB_type',1); % 1 or 2
stepsizeQN  = setOpts('quasiNewton_stepsize', .8*2/Lip );
S_L_S       = setOpts('quasiNewton_SLS', true );
displayTime = setOpts('displayTime',LARGE | HUGE );
SVDstyle    = setOpts('SVDstyle', 1*SMALL + 4*(~SMALL) ); % 1 is full SVD
% and even finer tuning (only matter if SVDstyle==4)
SVDwarmstart= setOpts('SVDwarmstart', true );
SVDnPower   = setOpts('SVDnPower', 1 + ~SVDwarmstart ); % number of power iteratiosn
SVDoffset   = setOpts('SVDoffset', 5 );
SVDopts = struct('SVDstyle', SVDstyle,'warmstart',SVDwarmstart,...
    'nPower',SVDnPower,'offset',SVDoffset );
if tau < 0
	Lagrangian = true;
	if sumProject || maxProject
		error('in Lagrangian mode (when tau<0 significies lambda=|tau|), turn off sum/maxProject');
	end
	lambda = abs(tau);
	tau = []; % help us track down bugs
else
	Lagrangian = false;
	if (sumProject && maxProject) || (~sumProject && ~maxProject), error('must choose either "sum" or "max" type projection'); end
end
if QUASINEWTON 
    if sumProject
        error('Can not run quasi-Newton mode when in "sum" formulation. Please change to "max"');
    elseif FISTA
        error('Can not run quasi-Newton with FISTA');
    elseif BB
        error('Can not run quasi-Newton with BB');
    end
end

projNuclear(); % remove any persistent variables
if maxProject
    project = @(L,S,varargin) projectMax(tau,lambda_S,SVDopts, L,S);
elseif sumProject
    project = @(L,S,varargin) projectSum(tau,lambda_S,L,S);
elseif Lagrangian
% 	project = @(L,S,varargin) proximityLS(lambda,lambda*lambda_S, L, S, varargin{:});
    project = @(L,S,varargin) projectMax(lambda,lambda*lambda_S,SVDopts, L,S, varargin{:});
end

L           = setOpts('L0',zeros(n1,n2) );
S           = setOpts('S0',zeros(n1,n2) );

% Check for extra options that were not processed
if ~isempty( fieldnames(opts ) )
    disp( 'warning, found extra guys in opts');
    disp( opts )
    error('Found unprocessed options in "opts"');
end

stepsize    = 1/Lip;
errHist     = zeros(maxIts, 2 + ~isempty(errFcn) );
Grad        = 0;
if FISTA || BB || QUASINEWTON
    L_old   = L;
    S_old   = S;
end
L_fista     = L;
S_fista     = S;
BREAK   = false;
kk      = 0; % counter for FISTA
timeRef = tic;
for k = 1:maxIts
    % Gradient in (L,S) (at the fista point) is (R,R) where...
    R   = A(L_fista + S_fista) - AY;
    Grad_old    = Grad;
    Grad        = At(R);
    
    if QUASINEWTON
%         stepsizeQN  = 1 - min(0.1, 1/k );
%         stepsizeQN = 1 - .3/sqrt(k);
        
        if S_L_S
            % we solve for S, update L, then re-update S
            % Exploits the fact that projection for S is faster
            dL      = L - L_old;
            S_old   = S;
            [~,S_temp]   = project( [], S - stepsizeQN*( Grad + dL ), stepsizeQN ); % take small step...
            
            dS      = S_temp - S_old;
            L_old   = L;
            [L,~,rnk]   = project( L - stepsizeQN*( Grad + dS ), [], stepsizeQN );
            
            dL      = L - L_old;
            [~,S]   = project( [], S - stepsizeQN*( Grad + dL ) , stepsizeQN);
        else
            % Gauss-Seidel update, starting with L, then S
            % Changing order seemed to not work as well
            dS      = S - S_old;
            L_old   = L;
            [L,~,rnk]   = project( L - stepsizeQN*( Grad + dS ), [] , stepsizeQN);
            dL      = L - L_old;
            S_old   = S;
            [~,S]   = project( [], S - stepsizeQN*( Grad + dL ) , stepsizeQN);
        end
%         stepsizeQN = (1+3*stepsizeQN)/4;
    else
        if BB && k > 1
            [stepsizeL, stepsizeS]  = compute_BB_stepsize( Grad, Grad_old, L, L_old, S, S_old, BB_split, BB_type);
            if isnan(stepsizeL) || isnan(stepsizeS)
                fprintf(2,'Warning: no BB stepsize possible since iterates have not changed!\n');
                [stepsizeL, stepsizeS]   = deal( stepsize );
            end
        else
            [stepsizeL, stepsizeS]   = deal( stepsize );
        end
        
        % Now compute proximity step
        if FISTA || BB
            L_old   = L;
            S_old   = S;
        end
        L           = L_fista - stepsizeL*Grad;
        S           = S_fista - stepsizeS*Grad;
        if any(isnan(L(:))) || any(isnan(S(:))), fprintf(2,'DEBUG!\n'); keyboard; end
        [L,S,rnk]   = project( L, S, stepsizeL, stepsizeS);
        if any(isnan(L(:))) || any(isnan(S(:))), fprintf(2,'DEBUG!\n'); keyboard; end
    end
    
    
    DO_RESTART = false;
    if FISTA
        if k>1 && restart > 0 && ~isinf(restart) && ~mod(kk,restart)
            kk = 0;
            DO_RESTART  = true;
        elseif restart==-Inf && kk > 5
            % In this case, we restart if the function has significantly increased
            if (errHist(k-1,2)-errHist(k-5,2)) > 1e-8*abs(errHist(k-5,2))
                DO_RESTART = true;
                kk = 0;
            end
        end
        L_fista = L + kk/(kk+3)*( L - L_old );
        S_fista = S + kk/(kk+3)*( S - S_old );
        kk      = kk + 1;
    else
        L_fista = L;
        S_fista = S;
    end

    
%     res          = norm(R(:)); % this is for L_fista, not L
    R            = A(L + S) - AY;
    % sometimes we already have R pre-computed, so if this turns out to 
    %   be a signficant computational cost we can make fancier code...
    res          = norm(R(:));
    errHist(k,1) = res;
    errHist(k,2) = 1/2*(res^2); % + lambda_L*objL + lambda_S*objS;
    if k > 1 && abs(diff(errHist(k-1:k,1)))/res < tol
        BREAK = true;
    end
    PRINT   = ~mod(k,printEvery) | BREAK | DO_RESTART;
    if PRINT
        fprintf('Iter %4d, residual %.2e, objective %.2e', k, res, errHist(k,2) -trueObj);
    end
    if ~isempty(errFcn)
        err     = errFcn(L,S);
        errHist(k,3)    = err;
        if PRINT, fprintf(', err %.2e', err ); end
    end
    if ~isempty(rnk) && PRINT, fprintf(', rank(L) %3d', rnk ); end
    if displayTime && PRINT
        tm = toc( timeRef ); timeRef = tic;
        fprintf(', time %.1f s', tm );
    end
    if DO_RESTART, fprintf(' [restarted FISTA]'); end
    if PRINT, fprintf('\n'); end
    if BREAK
        break;
    end
end
if BREAK, errHist = errHist(1:k,:); end

end

% subfunctions for projection
function [L,S,rnk] = projectMax( tau, lambda_S,SVDopts, L, S , stepsize, stepsizeS)
if nargin >= 6 && ~isempty(stepsize)
	% we compute proximity, not projection
	tauL 	= -abs( tau*stepsize );
    if nargin < 7 || isempty( stepsizeS ), stepsizeS = stepsize; end
	tauS 	= -abs( lambda_S*stepsizeS );
else
	tauL 	= abs(tau);
	tauS 	= abs(tau/lambda_S );
end

  % We project separately, so very easy
  if ~isempty(L)
%       projNuclear     = project_nuclear(tau); % in proxes/ folder
      [L,rnk]           = projNuclear(tauL, L,SVDopts);
      
      % for debugging
%       [L1,rnk]           = projNuclear(tauL, L,SVDopts);
%       [U,D,V] 	= svd(L,'econ');
%       d = max( 0, diag(D) - abs(tauL) );
%       rnk2 = nnz(d);
%       L2 = U*diag(d)*V';
%       if rnk ~= rnk2 || norm(L1-L2,'fro')/max(norm(L1,'fro'),norm(L2,'fro')) > 1e-3
%           fprintf(2,'Debugging!\n');
%           keyboard;
%       end
%       L=L1;
      
  else
      rnk = [];
  end
  if ~isempty(S)
	  if tauS > 0
        projL1          = project_l1(tauS);
        S               = projL1(      S);
	  else
		% simple prox
		S = sign(S).*max(0, abs(S) - abs(tauS));
	  end
  end
end
function [X,rEst] = projNuclear( tau, X, SVDopts )
% Input must be a matrix, not a vector
persistent oldRank Vold iteration
if nargin==0, oldRank=[]; Vold = []; iteration = 0; return; end
if isempty(oldRank), rEst = 10;
else, rEst = oldRank + 2;
end
if isempty(iteration), iteration = 0; end
iteration   = iteration + 1;
[n1,n2]     = size(X);
minN       = min( [n1,n2] ); % could set smaller to make nonconvex
% For the first few iterations, we constrain rankMax
switch iteration
    case 1
        rankMax     = round(minN/4);
    case 2
        rankMax     = round(minN/2);
    otherwise
        rankMax     = minN;
end

style = SVDopts.SVDstyle;
if tau==0, X=0*X; return; end

switch style
    case 1
        % full SVD
        [U,S,V] = svd(X,'econ');
        s   = diag(S);
        if tau < 0 % we do prox
            s = max( 0, s - abs(tau) );
        else
            s = project_simplex( tau, s );
        end
        tt      = s > 0;
        rEst    = nnz(tt);
        U       = U(:,tt);
        S       = diag(s(tt));
        V       = V(:,tt);
    case {2,3,4}
        % 2: use Matlab's sparse SVD
        % 3: use PROPACK
        % 4: use Joel Tropp's randomized SVD
        if style==2
            opts = struct('tol',1e-4);
            if rankMax==1, opts.tol = min(opts.tol,1e-6); end % important!
            svdFcn = @(X,rEst)svds(X,rEst,'L',opts); 
        elseif style==3
            opts = struct('tol',1e-4,'eta',eps);
            opts.delta = 10*opts.eta;
            % set eta to eps, but not 0 otherwise reorth is very slow
            if rankMax==1, opts.tol = min(opts.tol,1e-6); end % important!
            svdFcn = @(X,rEst)lansvd(X,rEst,'L',opts); 
        elseif style == 4
            opts = [];
            if isfield(SVDopts,'nPower') && ~isempty( SVDopts.nPower )
                nPower = SVDopts.nPower;
            else
                nPower = 2;
            end
            if isfield(SVDopts,'offset') && ~isempty( SVDopts.offset )
                offset = SVDopts.offset;
            else
                offset = 5;
            end
            if isfield( SVDopts, 'warmstart' ) && SVDopts.warmstart==true ...
                    && ~isempty(Vold)
                opts = struct( 'warmStart', Vold );
            end
            ell     = @(r) min([r+offset,n1,n2]); % number of samples to take
            svdFcn = @(X,rEst)randomizedSVD(X,rEst,ell(rEst),nPower,[],opts );
        end
        
        ok  = false;
        while ~ok
            rEst    = min( [rEst,rankMax] );
            [U,S,V] = svdFcn(X,rEst);
            s       = diag(S);
            if tau < 0
                % we are doing prox
                lambda = abs(tau);
            else
                lambda  = findTau(s,tau);
            end
            ok      = ( min(s) < lambda ) || (rEst == rankMax);
            if ok, break; end
            rEst    = 2*rEst;
        end
        rEst = min( length(find(s>lambda)), rankMax );
        S   = diag( s(1:rEst) - lambda );
        U   = U(:,1:rEst);
        V   = V(:,1:rEst);
    otherwise
        error('bad value for SVDstyle');
end
if isempty(U)
    X = 0*X;
else
    X = U*S*V';
end
oldRank = size(U,2);
if isfield( SVDopts, 'warmstart' ) && SVDopts.warmstart==true
    Vold = V;
end
end

function x = project_simplex( q, x )
% projects onto the constraints sum(x)=q and x >= 0
% Update: projects onto sum(x) <= q and x >= 0

x     = x.*( x > 0 ); % March 11 '14
% March 11, fixing bug: we want to project onto the volume, not
%   the actual simplex (surface)
if sum(x) <= q, return; end
if q==0, x = 0*x; return; end

s     = sort( x, 'descend' );
if q < eps(s(1))
    % eps(x) is the distance from abs(x) to the next larger
    %   floating point number, i.e., in floating point arithmetic,
    %   x + eps(x)/2 = x
    % eps(1) is about 2.2e-16
    
    error('Input is scaled so large compared to q that accurate computations are difficult');
    % since then cs(1) = s(1) - q is  not even guaranteed to be
    % smaller than q !
end
cs    = ( cumsum(s) - q ) ./ ( 1 : numel(s) )';
ndx   = nnz( s > cs );
x     = max( x - cs(ndx), 0 );
end

function tau = findTau( s, lambda )
% Returns the shrinkage value necessary to shrink the vector s
%   so that it is in the lambda scaled simplex
%   Usually, s is the diagonal part of an SVD
if all(s==0)||lambda==0, tau=0; return; end
if numel(s)>length(s), error('s should be a vector, not a matrix'); end
if ~issorted(flipud(s)), s = sort(s, 'descend'); end
if any( s < 0 ), error('s should be non-negative'); end

% project onto the simplex of radius lambda (not tau)
% and use this to find "tau" (the shrinkage amount)
% If we know s_i > tau for i = 1, ..., k, then
%   tau = ( sum(s(1:k)) - lambda )/k
% But we don't know k, so find it:
cs  = (cumsum(s) - abs(lambda) )./(1:length(s))';
ndx = nnz( s > cs ); % >= 1 as long as lambda > 0
tau = max(0,cs(ndx));
% We want to make sure we project onto sum(sigma) <= lambda
%   and not sum(sigma) == lambda, so do not allow negative tau

end


function [L,S,rnk] = projectSum( tau, lambda_S, L, S )
  [m,n]           = size(L);
  [U,Sigma,V]     = svd(L,'econ');
  s       = diag(Sigma);
  wts     = [ ones(length(s),1); lambda_S*ones(m*n,1) ];
  proj    = project_l1(tau,wts);
  sS      = proj( [s;vec(S)] );
  sProj   = sS(1:length(s));
  S       = reshape( sS(length(s)+1:end), m, n );
  L       = U*diag(sProj)*V';
  rnk     = nnz( sProj );
end


function [stepsizeL, stepsizeS]  = compute_BB_stepsize( ...
    Grad, Grad_old, L, L_old, S, S_old, BB_split, BB_type)

  if ~BB_split
      % we take a Barzilai-Borwein stepsize in the full variable
      yk  = Grad(:) - Grad_old(:);
      yk  = [yk;yk]; % to account for both variables
      sk  = [L(:) - L_old(:); S(:) - S_old(:) ];
      if BB_type == 1
          % Default. The bigger stepsize
          stepsize    = norm(sk)^2/(sk'*yk);
      elseif BB_type == 2
          stepsize    = sk'*yk/(norm(yk)^2);
      end
      [stepsizeL, stepsizeS] = deal( stepsize );
      
  elseif BB_split
      % treat L and S variables separately
      % Doesn't seem to work well.
      yk  = Grad(:) - Grad_old(:);
      skL  = L(:) - L_old(:);
      skS  = S(:) - S_old(:);
      if BB_type == 1
          % Default. The bigger stepsize
          stepsizeL   = norm(skL)^2/(skL'*yk);
          stepsizeS   = norm(skS)^2/(skS'*yk);
      elseif BB_type == 2
          stepsizeL   = skL'*yk/(norm(yk)^2);
          stepsizeS   = skS'*yk/(norm(yk)^2);
      end
  end
end

function [L,S,rnk] = proximityLS( lambdaL, lambdaS, L, S, stepsizeL, stepsizeS )
if nargin < 6 || isempty( stepsizeS )
    stepsizeS = stepsizeL;
end
rnk = [];
if ~isempty(L)
	[U,D,V] 	= svd(L,'econ');
	d = max( 0, diag(D) - (lambdaL*stepsizeL) );
	rnk = nnz(d);
	L = U*diag(d)*V';
end
if ~isempty(S)
	S = sign(S).*max( 0, abs(S) - (lambdaS*stepsizeS) );
end
end






function op = project_l1( q , d)
%PROJECT_L1   Projection onto the scaled 1-norm ball.
%    OP = PROJECT_L1( Q ) returns an operator implementing the 
%    indicator function for the 1-norm ball of radius q,
%    { X | norm( X, 1 ) <= q }. Q is optional; if omitted,
%    Q=1 is assumed. But if Q is supplied, it must be a positive
%    real scalar.
%
%    OP = PROJECT_L1( Q, D ) uses a scaled 1-norm ball of radius q,
%    { X | norm( D.*X, 1 ) <= 1 }. D should be the same size as X
%    and non-negative (some zero entries are OK).

% Note: theoretically, this can be done in O(n)
%   but in practice, worst-case O(n) median sorts are slow
%   and instead average-case O(n) median sorts are used.
%   But in matlab, the median is no faster than sort
%   (the sort is probably quicksort, O(n log n) expected, with
%    good constants, but O(n^2) worst-case).
%   So, we use the naive implementation with the sort, since
%   that is, in practice, the fastest.

if nargin == 0,
	q = 1;
elseif ~isnumeric( q ) || ~isreal( q ) || numel( q ) ~= 1 || q <= 0,
	error( 'Argument must be positive.' );
end
if nargin < 2 || isempty(d) || numel(d)==1
    if nargin>=2 && ~isempty(d)
        % d is a scalar, so norm( d*x ) <= q is same as norm(x)<=q/d
        if d==0
            error('If d==0 in proj_l1, the set is just {0}, so use proj_0');
        elseif d < 0
            error('Require d >= 0');
        end
        q = q/d;
    end
    op = @(varargin)proj_l1_q(q, varargin{:} );
else
    if any(d<0)
        error('All entries of d must be non-negative');
    end
    op = @(varargin)proj_l1_q_d(q, d, varargin{:} );
end

% This is modified from TFOCS, Nov 26 2013, Stephen Becker
% Note: removing "v" (value) output from TFOCS code
%   Also removing extraneous inputs
function x = proj_l1_q( q, x, varargin )
    myReshape   = @(x) x; % do nothing
    if size(x,2) > 1
        if ndims(x) > 2, error('You must modify this code to deal with tensors'); end
        myReshape     = @(y) reshape( y, size(x,1), size(x,2) );
        x   = x(:); % make it into a vector
    end
    s      = sort(abs(nonzeros(x)),'descend');
    cs     = cumsum(s);
    % ndx    = find( cs - (1:numel(s))' .* [ s(2:end) ; 0 ] >= q, 1 );
    ndx    = find( cs - (1:numel(s))' .* [ s(2:end) ; 0 ] >= q+2*eps(q), 1 ); % For stability
    if ~isempty( ndx )
        thresh = ( cs(ndx) - q ) / ndx;
        x      = x .* ( 1 - thresh ./ max( abs(x), thresh ) ); % May divide very small numbers
    end
    x   = myReshape(x);
end

% Allows scaling. Added Feb 21 2014
function x = proj_l1_q_d( q, d,  x, varargin )
    myReshape   = @(x) x; % do nothing
    if size(x,2) > 1
        if ndims(x) > 2, error('You must modify this code to deal with tensors'); end
        myReshape     = @(y) reshape( y, size(x,1), size(x,2) );
        x   = x(:); % make it into a vector
    end
    [goodInd,j,xOverD] = find( x./ d );
    [lambdas,srt]      = sort(abs(xOverD),'descend');
    s   = abs(x(goodInd).*d(goodInd));  
    s   = s(srt);
    dd  = d(goodInd).^2;
    dd  = dd(srt);
    cs  = cumsum(s);
    cd  = cumsum(dd);
    ndx    = find( cs - lambdas.*cd >= q+2*eps(q), 1, 'first');
    if ~isempty( ndx )
        ndx     = ndx - 1;
        lambda  = ( cs(ndx) - q )/cd(ndx);
        x       = sign(x).*max( 0, abs(x) - lambda*d );
    end
    x   = myReshape(x);
end


end
