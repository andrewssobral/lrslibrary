function [T, max_mem, Uinit] = tucker_me(X, R, esz, opts)
%TUCKER_ME Memory-efficient Tucker higher-order orthogonal iteration.
%
%   T = TUCKER_ME(X,R,ESZ) computes the best rank(R1,R2,..,Rn)
%   approximation of tensor X, according to the specified dimensions
%   in vector R. ESZ specifies the number of dimensions that are to be
%   handled elementwise. The input X is sptensor. The result
%   returned in T is a ttensor. 
%
%   T = TUCKER_ME(X,R,ESZ,OPTS) specify options:
%   OPTS.tol: Tolerance on difference in fit {1.0e-4}
%   OPTS.maxiters: Maximum number of iterations {50}
%   OPTS.dimorder: Order to loop through dimensions {1:ndims(A)}
%   OPTS.init: Initial guess [{'random'}|'eigs'|cell array]
%
%   [T,MEM,U0] = TUCKER_ME(...) also returns the memory estimate and
%   the initial guess U0.  
%
%   See also TENSOR_TOOLBOX, TUCKER_ALS, TTM_ME, TTM_ME_PARTITION,
%   TTM_ME_MEM.
%
%   Code by Tamara Kolda and Jimeng Sun, 2008. Adapted with permission from
%   the MATLAB Tensor Toolbox function algorithms/tucker_als.
%
%   Based on the paper:
%   T. G. Kolda and J. Sun. Scalable Tensor Decompositions for Multi-aspect
%   Data Mining. In: ICDM 2008: Proceedings of the 8th IEEE International
%   Conference on Data Mining, December 2008.  

% $Id: tucker_me.m,v 1.3 2010/03/19 22:55:46 tgkolda Exp $


%% Fill in optional variable
if ~exist('opts','var')
    opts = struct;
end

%% Extract number of dimensions and norm of X.
N = ndims(X);
normX = norm(X);

%% Set algorithm parameters from input or by using defaults
fitchangetol = setparam(opts,'tol',1e-4);
maxiters = setparam(opts,'maxiters',50);
dimorder = setparam(opts,'dimorder',1:N);
init = setparam(opts,'init','random');

if numel(R) == 1
    R = R * ones(N,1);
end

%% Error checking 
% Error checking on maxiters
if maxiters < 0
    error('OPTS.maxiters must be positive');
end

% Error checking on dimorder
if ~isequal(1:N,sort(dimorder))
    error('OPTS.dimorder must include all elements from 1 to ndims(X)');
end

%% Set up and error checking on initial guess for U.
if iscell(init)
    Uinit = init;
    if numel(Uinit) ~= N
        error('OPTS.init does not have %d cells',N);
    end
    for n = dimorder(2:end);
        if ~isequal(size(Uinit{n}),[size(X,n) R(n)])
            error('OPTS.init{%d} is the wrong size',n);
        end
    end
else
    % Observe that we don't need to calculate an initial guess for the
    % first index in dimorder because that will be solved for in the first
    % inner iteration.
    if strcmp(init,'random')
        Uinit = cell(N,1);
        for n = dimorder(2:end)
            Uinit{n} = rand(size(X,n),R(n));
        end
    elseif strcmp(init,'nvecs') || strcmp(init,'eigs') 
        % Compute an orthonormal basis for the dominant
        % Rn-dimensional left singular subspace of
        % X_(n) (1 <= n <= N).
        Uinit = cell(N,1);
        for n = dimorder(2:end)
            fprintf('  Computing %d leading e-vectors for factor %d.\n', ...
                    R(n),n);
            Uinit{n} = nvecs(X,n,R(n));
        end
    else
        error('The selected initialization method is not supported');
    end
end

%% Set up for iterations - initializing U, fit, and max_mem
U = Uinit;
fit = 0;
max_mem = 0;

%% Main Loop: Iterate until convergence
fprintf('\nAlternating Least-Squares:\n');
for iter = 1:maxiters
    fitold = fit;

    % Iterate over all N modes of the tensor
    for n = dimorder(1:end)
        % Estimate the best partition of edims and sdims based on n
        [ed, sd] = ttm_me_partition(U, esz, n);
    	Utilde = ttm_me(X, U, ed, sd, 't');
        % Check memory
        mem = ttm_me_mem(X, U, ed, sd, 't');
        max_mem = max(mem, max_mem);
        % Maximize norm(Utilde x_n W') wrt W and
        % keeping orthonormality of W
        U{n} = nvecs(Utilde,n,R(n));
    end

    % Assemble the current approximation
    core = ttm(Utilde, U, n,  't');  

    % Compute fit
    normresidual = sqrt( normX^2  - norm(core)^2);
    fit = 1 - (normresidual / normX); % fraction explained by model
    fitchange = abs(fitold - fit);

    fprintf(' Iter %2d: fit = %e fitdelta = %7.1e\n', iter, fit, fitchange);

    % Check for convergence
    if (iter > 1) && (fitchange < fitchangetol)
        break;
    end

end

%% Compute the final result %commented by Jimeng Sun, seems redundant
max_mem = max(mem, max_mem);    

%% Assemble the resulting tensor
T = ttensor(core, U);
end

%%
function x = setparam(opts,name,default)
if isfield(opts,name);
    x = opts.(name);
else
    x = default;
end
end
%%
 