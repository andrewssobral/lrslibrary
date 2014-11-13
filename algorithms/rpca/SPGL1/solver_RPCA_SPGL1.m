function [L,S,errHist] = solver_RPCA_SPGL1(AY,lambda_S, epsilon, A_cell, opts)
% [L,S,errHist] = solver_RPCA_SPGL1(Y,lambda_S, epsilon, A_cell, opts)
% Solves the problem
%   minimize_{L,S} f(L,S)
%   subject to 
%       || A( L+S ) - Y ||_F <= epsilon
%   if opts.sum = true
%       (1) f(L,S) = ||L||_* + lambda_S ||S||_1 
%   if opts.max = true
%       (2) f(L,S) = max( ||L||_* , lambda_S ||S||_1 ) 
%
%   By default, A is the identity, 
%   or if A_cell is provided, where A_cell = {A, At}
%   (A is a function handle, At is a function handle to the transpose of A)
%
%   This uses the trick from SPGL1 to call several subproblems
%
%   errHist(1,:) is a record of the residual
%   errHist(2,:) is a record of the full objective (that is, .5*resid^2 )
%   errHist(3,:) is the output of opts.errFcn if provided
%
% opts is a structure with options:
%   opts.tau0       starting point for SPGL1 search (default: 10)
%   opts.SPGL1_tol  tolerance for 1D Newton search (default: 1e-2)
%   opts.SPGL1_maxIts   maximum number of iterations for 1D search (default: 10)
%
% And these options are the same as in solver_RPCA_constrained.m
%   opts.sum, opts.max  (as described above)
%   opts.L0         initial guess for L (default is 0)
%   opts.S0         initial guess for S (default is 0)
%   opts.tol        sets stopping tolerance (default is 1e-6)
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
% See also solver_RPCA_Lagrangian.m, solver_RPCA_constrained.m


% todo: allow S >= 0 constraints, since this is easy

error(nargchk(3,5,nargin,'struct'));
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
else
    A   = A_cell{1};
    At  = A_cell{2};
    % Y could be either Y or A(Y)
    if size(AY,2) > 1
        % AY is a vector, so it is probably Y and not AY
        disp('Changing Y to A(Y)');
        AY = A(AY);
    end
end

tauInitial      = setOpts('tau0', 1e1 );


finalTol        = setOpts('tol',1e-6);
sumProject      = setOpts('sum', false );
maxProject      = setOpts('max', false );
if (sumProject && maxProject) || (~sumProject && ~maxProject), error('must choose either "sum" or "max" type projection'); end
opts.max        = maxProject; % since setOpts removes them
opts.sum        = sumProject;
rNormOld        = Inf;
tau             = tauInitial;
SPGL1_tol       = setOpts('SPGL1_tol',1e-2);
SPGL1_maxIts    = setOpts('SPGL1_maxIts', 10 );
errHist         = [];
for nSteps = 1:SPGL1_maxIts
    opts.tol    = max(finalTol, finalTol*10^((4/nSteps)-1) ); % slowly reduce tolerance
    %         opts.tol    = finalTol;
    fprintf('\n==Running SPGL1 Newton iteration with tau=%.2f\n\n', tau);
    [L,S,errHistTemp] = solver_RPCA_constrained(AY,lambda_S, tau, A_cell, opts);
    % Find a zero of the function phi(tau) = norm(x_tau) - epsilon
    
    errHist = [errHist; errHistTemp];
    rNorm       = errHist(end,1); % norm(residual)
    if abs( epsilon - rNorm ) < SPGL1_tol*epsilon
        disp('Reached end of SPGL1 iterations: converged to right residual');
        break;
    end
    if abs( rNormOld - rNorm ) < .1*SPGL1_tol
        disp('Reached end of SPGL1 iterations: converged');
        break;
    end
    % for nuclear norm matrix completion,
    %   normG = spectralNorm( gradient ) and gradient is just residual
    G   = At(A(L+S-AY)); % it's really two copies, [G;G]
    % if we have a linear operator, above needs to be modified...
    if sumProject
        normG = max(norm(G), (1/lambda_S)*norm(G(:), inf));
    elseif maxProject
        normG = norm(G) +  (1/lambda_S)*norm(G(:), inf);
    end
    
    % otherwise, take another Newton step
    phiPrime    = -normG/rNorm; % derivative of phi
    ALPHA       = .99; % amount of Newton step to take
    tauOld      = tau;
    tau         = tau + ALPHA*( epsilon - rNorm )/phiPrime;
    tau         = min( max(tau,tauOld/10), 1e10 );
    % and warm-start the next iteration:
    opts.S0     = S;
    opts.L0     = L;
end
% And if we requested very high accuracy, run a bit more on the final
% iteration...
if finalTol < 1e-6
    opts.tol    = finalTol/10;
    opts.S0     = S;
    opts.L0     = L;
    opts.SVDnPower  = 3;
    [L,S,errHistTemp] = solver_RPCA_constrained(AY,lambda_S, tau, A_cell, opts);
    errHist = [errHist; errHistTemp];
end


end % end of main function