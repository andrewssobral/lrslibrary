function [M,Minit,output] = cp_apr(X, R, varargin)
%CP_APR Compute nonnegative CP with alternating Poisson regression.
%
%   M = CP_APR(X, R) computes an estimate of the best rank-R
%   CP model of a tensor X using an alternating Poisson regression.
%   The input X can be a tensor, sptensor, ktensor, or ttensor. The
%   result P is a ktensor.
%
%   M = CP_APR(X, R, 'param', value, ...) specifies optional parameters and
%   values. Valid parameters and their default values are:
%      'tol' - Tolerance on the inner KKT violation {1.0e-4}
%      'maxiters' - Maximum number of iterations {1000}
%      'maxinneriters' = Maximum number of inner iterations {10}
%      'init' - Initial guess [{'random'}|ktensor]
%      'epsilon' - parameter to avoid divide by zero {100*eps}
%      'kappatol' - tolerance on complementary slackness {100*eps}
%      'kappa' - offset to fix complementary slackness {10*eps}
%      'printitn' - Print every n outer iterations; 0 for no printing {1}
%      'printinneritn' - Print every n inner iterations {0}
%
%   [M,M0] = CP_APR(...) also returns the initial guess.
%
%   [M,M0,out] = CP_APR(...) also returns additional output.
%      out.kktViolations - maximum kkt violation per iteration
%      out.nInnerIters   - number of inner iterations per iteration
%      out.nViolations   - number of factor matrices needing complementary
%                          slackness adjustment per iteration
%      out.nTotalIters   - total number of inner iterations
%
%   REFERENCE: E. C. Chi and T. G. Kolda. On Tensors, Sparsity, and
%   Nonnegative Factorizations, arXiv:1112.2414 [math.NA], December 2011,
%   URL: http://arxiv.org/abs/1112.2414. Submitted for publication.
%
%   See also CP_ALS, KTENSOR, TENSOR, SPTENSOR.
%
%MATLAB Tensor Toolbox.
%Copyright 2012, Sandia Corporation.

% This is the MATLAB Tensor Toolbox by T. Kolda, B. Bader, and others.
% http://www.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2012) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in the file LICENSE.txt

%% Extract dimensions of X and number of dimensions of X.
N = ndims(X);

%% Set algorithm parameters from input or by using defaults.
params = inputParser;
params.addParamValue('epsilon',1e-10,@isscalar);
params.addParamValue('tol',1e-4,@isscalar);
params.addParamValue('maxiters',1000,@(x) isscalar(x) & x > 0);
params.addParamValue('init','random',@(x) (isa(x,'ktensor') || ismember(x,{'random'})));
params.addParamValue('printitn',1,@isscalar);
params.addParamValue('kappa',1e-2,@isscalar);
params.addParamValue('kappatol',1e-10,@isscalar);
params.addParamValue('maxinneriters',10,@isscalar);
params.addParamValue('printinneritn',0,@isscalar);
params.parse(varargin{:});


%% Copy from params object.
epsilon = params.Results.epsilon;
tol = params.Results.tol;
maxOuterIters = params.Results.maxiters;
Minit = params.Results.init;
kappa = params.Results.kappa;
kappaTol = params.Results.kappatol;
maxInnerIters = params.Results.maxinneriters;
printOuterItn = params.Results.printitn;
printInnerItn = params.Results.printinneritn;
kktViolations = -ones(maxOuterIters,1);
nInnerIters = zeros(maxOuterIters,1);

%% Set up and error checking on initial guess for U.
if isa(Minit,'ktensor')
    if ndims(Minit) ~= N
        error('Initial guess does not have the right number of dimensions');
    end
    
    if ncomponents(Minit) ~= R
        error('Initial guess does not have the right number of components');
    end
    
    for n = 1:N
        if size(Minit,n) ~= size(X,n)
            error('Dimension %d of the initial guess is the wrong size',n);
        end
    end
elseif strcmp(Minit,'random')
    F = cell(N,1);
    for n = 1:N
        F{n} = rand(size(X,n),R);
    end
    Minit = ktensor(F);
else
    error('The selected initialization method is not supported');
end


%% Set up for iterations - initializing M and Phi.
M = normalize(Minit,[],1);
Phi = cell(N,1);
kktModeViolations = zeros(N,1);

if printOuterItn > 0
  fprintf('\nCP_APR:\n');
end

nViolations = zeros(maxOuterIters,1);

%% Main Loop: Iterate until convergence.
for iter = 1:maxOuterIters
    
    isConverged = true;   
    for n = 1:N

        % Make adjustments to entries of M{n} that are violating
        % complementary slackness conditions.
        if (iter > 1)
            V = (Phi{n} > 1) & (M{n} < kappaTol);
            if any(V(:))           
                nViolations(iter) = nViolations(iter) + 1;
                M{n}(V>0) = M{n}(V>0) + kappa;
            end
        end         

        % Shift the weight from lambda to mode n
        M = redistribute(M,n);
        
        % Calculate product of all matrices but the n-th
        % (In sparse case, only calcuates entries corresponding to nonzeros in X.)
        Pi = calculatePi(X, M, R, n, N);        
        
        % Do the multiplicative updates
        for i = 1:maxInnerIters

            % Count the inner iterations
            nInnerIters(iter) = nInnerIters(iter) + 1;
                                  
            % Calculate matrix for multiplicative update
            Phi{n} = calculatePhi(X, M, R, n, Pi, epsilon);
            
            % Check for convergence
            kktModeViolations(n) = max(abs(vec(min(M.U{n},1-Phi{n}))));
            if (kktModeViolations(n) < tol)
                break;
            else
                isConverged = false;
            end                      
            
            % Do the multiplicative update
            M{n} = M{n} .* Phi{n};
            
            % Print status
             if mod(i, printInnerItn)==0
                 fprintf('    Mode = %1d, Inner Iter = %2d, KKT violation = %.6e\n', n, i, kktModeViolations(n));
             end
        end
        
        % Shift weight from mode n back to lambda
        M = normalize(M,[],1,n);
        
    end

    kktViolations(iter) = max(kktModeViolations);    

    if (mod(iter,printOuterItn)==0)
        fprintf(' Iter %4d: Inner Its = %2d KKT violation = %.6e, nViolations = %2d\n', ...
        iter, nInnerIters(iter), kktViolations(iter), nViolations(iter));            
    end
    
    % Check for convergence
    if (isConverged)
        break;
    end    
end

%% Clean up final result
M = normalize(M,'sort',1);

if printOuterItn>0
    normX = norm(X);   
    normresidual = sqrt( normX^2 + norm(M)^2 - 2 * innerprod(X,M) );
    fit = 1 - (normresidual / normX); %fraction explained by model
    fprintf('===========================================\n');
    fprintf(' Final log-likelihood = %e \n', tt_loglikelihood(X,M));
    fprintf(' Final least squares fit = %e \n', fit);
    fprintf(' Final KKT violation = %7.7e\n', kktViolations(iter));
    fprintf(' Total inner iterations = %d\n', sum(nInnerIters));
end

output = struct;
output.params = params.Results;
output.kktViolations = kktViolations(1:iter);
output.nInnerIters = nInnerIters(1:iter);
output.nViolations = nViolations(1:iter);
output.nTotalIters = sum(nInnerIters);


end

function Pi = calculatePi(X, M, R, n, N)

if (isa(X,'sptensor'))
    Pi = ones(nnz(X), R);
    for nn = [1:n-1,n+1:N]
        Pi = M{nn}(X.subs(:,nn),:).*Pi;
    end
else
    U = M.U;
    Pi = khatrirao(U{[1:n-1,n+1:N]},'r');
end

end

function Phi = calculatePhi(X, M, R, n, Pi, epsilon)

if (isa(X,'sptensor'))
    Phi = -ones(size(X,n),R);
    xsubs = X.subs(:,n);
    v = sum(M.U{n}(xsubs,:).*Pi,2);
    wvals = X.vals ./ max(v, epsilon);
    for r = 1:R
        Yr = accumarray(xsubs, wvals .* Pi(:,r), [size(X,n) 1]);
        Phi(:,r) = Yr;
    end    
else
    Xn = double(tenmat(X,n));
    V = M.U{n}*Pi';
    W = Xn ./ max(V, epsilon);
    Y = W * Pi;
    Phi = Y;
end

end

function y = vec(x)
y = x(:);
end
