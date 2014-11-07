function [T,Uinit] = tucker_als(X,R,varargin)
%TUCKER_ALS Higher-order orthogonal iteration.
%
%   T = TUCKER_ALS(X,R) computes the best rank-(R1,R2,..,Rn)
%   approximation of tensor X, according to the specified dimensions
%   in vector R.  The input X can be a tensor, sptensor, ktensor, or
%   ttensor.  The result returned in T is a ttensor.
%
%   T = TUCKER_ALS(X,R,'param',value,...) specifies optional parameters and
%   values. Valid parameters and their default values are:
%      'tol' - Tolerance on difference in fit {1.0e-4}
%      'maxiters' - Maximum number of iterations {50}
%      'dimorder' - Order to loop through dimensions {1:ndims(A)}
%      'init' - Initial guess [{'random'}|'nvecs'|cell array]
%      'printitn' - Print fit every n iterations {1}
%
%   [T,U0] = TUCKER_ALS(...) also returns the initial guess.
%
%   Examples:
%   X = sptenrand([5 4 3], 10);
%   T = tucker_als(X,2);        %<-- best rank(2,2,2) approximation 
%   T = tucker_als(X,[2 2 1]);  %<-- best rank(2,2,1) approximation 
%   T = tucker_als(X,2,'dimorder',[3 2 1]);
%   T = tucker_als(X,2,'dimorder',[3 2 1],'init','nvecs');
%   U0 = {rand(5,2),rand(4,2),[]}; %<-- Initial guess for factors of T
%   T = tucker_als(X,2,'dimorder',[3 2 1],'init',U0);
%
%   See also TTENSOR, TENSOR, SPTENSOR, KTENSOR.
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


% Extract number of dimensions and norm of X.
N = ndims(X);
normX = norm(X);

%% Set algorithm parameters from input or by using defaults
params = inputParser;
params.addParamValue('tol',1e-4,@isscalar);
params.addParamValue('maxiters',50,@(x) isscalar(x) & x > 0);
params.addParamValue('dimorder',1:N,@(x) isequal(sort(x),1:N));
params.addParamValue('init', 'random', @(x) (iscell(x) || ismember(x,{'random','nvecs','eigs'})));
params.addParamValue('printitn',1,@isscalar);
params.parse(varargin{:});

%% Copy from params object
fitchangetol = params.Results.tol;
maxiters = params.Results.maxiters;
dimorder = params.Results.dimorder;
init = params.Results.init;
printitn = params.Results.printitn;

if numel(R) == 1
    R = R * ones(N,1);
end
U = cell(N,1);

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

%% Set up for iterations - initializing U and the fit.
U = Uinit;
fit = 0;

if printitn > 0
    fprintf('\nTucker Alternating Least-Squares:\n');
end

%% Main Loop: Iterate until convergence
for iter = 1:maxiters

    fitold = fit;

    % Iterate over all N modes of the tensor
    for n = dimorder(1:end)
        Utilde = ttm(X, U, -n, 't');
        
        % Maximize norm(Utilde x_n W') wrt W and
        % keeping orthonormality of W
        U{n} = nvecs(Utilde,n,R(n));
    end

    % Assemble the current approximation
    core = ttm(Utilde, U, n, 't');

    % Compute fit
    normresidual = sqrt( normX^2 - norm(core)^2 );
    fit = 1 - (normresidual / normX); %fraction explained by model
    fitchange = abs(fitold - fit);

    if mod(iter,printitn)==0
        fprintf(' Iter %2d: fit = %e fitdelta = %7.1e\n', iter, fit, fitchange);
    end

    % Check for convergence
    if (iter > 1) && (fitchange < fitchangetol)
        break;
    end

end

T = ttensor(core, U);

end


