function [P,Uinit] = cp_nmu(X,R,opts)
%CP_NMU Compute nonnegative CP with multiplicative updates.
%
%   P = CP_NMU(X,R) computes an estimate of the best rank-R PARAFAC
%   model of a tensor X with nonnegative constraints on the factors.
%   This version uses the Lee & Seung multiplicative updates from
%   their NMF algorithm.  The input X can be a tensor, sptensor,
%   ktensor, or ttensor. The result P is a ktensor.
%
%   P = CP_NMU(X,R,OPTS) specify options:
%   OPTS.tol: Tolerance on difference in fit {1.0e-4}
%   OPTS.maxiters: Maximum number of iterations {50}
%   OPTS.dimorder: Order to loop through dimensions {1:ndims(A)}
%   OPTS.init: Initial guess [{'random'}|'nvecs'|cell array]
%   OPTS.printitn: Print fit every n iterations {1}
%
%   [P,U0] = CP_NMU(...) also returns the initial guess.
%
%   Examples:
%   X = sptenrand([5 4 3], 10);
%   P = cp_nmu(X,2);
%   P = cp_nmu(X,2,struct('dimorder',[3 2 1]));
%   P = cp_nmu(X,2,struct('dimorder',[3 2 1],'init','nvecs'));
%   U0 = {rand(5,2),rand(4,2),[]}; %<-- Initial guess for factors of P
%   P = cp_nmu(X,2,struct('dimorder',[3 2 1],'init',{U0}));
%
%   See also KTENSOR, TENSOR, SPTENSOR, TTENSOR.
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


%% Fill in optional variable
if ~exist('opts','var')
    opts = struct;
end

%% Extract number of dimensions and norm of X.
N = ndims(X);
normX = norm(X);

%% Set algorithm parameters from input or by using defaults
fitchangetol = setparam(opts,'tol',1e-4);
maxiters = setparam(opts,'maxiters',500);
dimorder = setparam(opts,'dimorder',1:N);
init = setparam(opts,'init','random');
printitn = setparam(opts,'printitn',1);
epsilon = 1e-12;  % Small number to protect against round-off error

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
    for n = dimorder(1:end);
        if ~isequal(size(Uinit{n}),[size(X,n) R])
            error('OPTS.init{%d} is the wrong size',n);
        end
    end
else
    if strcmp(init,'random')
        Uinit = cell(N,1);
        for n = dimorder(1:end)
            Uinit{n} = rand(size(X,n),R) + 0.1;
        end
    elseif strcmp(init,'nvecs') || strcmp(init,'eigs') 
        Uinit = cell(N,1);
        for n = dimorder(1:end)
            k = min(R,size(X,n)-2);
            fprintf('  Computing %d leading e-vectors for factor %d.\n',k,n);
            Uinit{n} = abs(nvecs(X,n,k));
            if (k < R)
              Uinit{n} = [Uinit{n} rand(size(X,n),R-k)]; 
            end
        end
    else
        error('The selected initialization method is not supported');
    end
end

%% Set up for iterations - initializing U and the fit.
U = Uinit;
fit = 0;

if printitn>0
  fprintf('\nNonnegative PARAFAC:\n');
end

%% Main Loop: Iterate until convergence
for iter = 1:maxiters

    fitold = fit;

    % Iterate over all N modes of the tensor
    for n = dimorder(1:end)

        % Compute the matrix of coefficients for linear system
        Y = ones(R,R);
        for i = [1:n-1,n+1:N]
            Y = Y .* (U{i}'*U{i});
        end
        Y = U{n} * Y;

        % Initialize matrix of unknowns
        Unew = U{n};

        % Calculate Unew = X_(n) * khatrirao(all U except n, 'r').
        tmp = mttkrp(X,U,n) + epsilon;

        % Update unknowns
        Unew = Unew .* tmp;
        Unew = Unew ./ (Y + epsilon);

        U{n} = Unew;
    end

    P = ktensor(U);
    normresidual = sqrt( normX^2 + norm(P)^2 - 2 * innerprod(X,P) );
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

%% Clean up final result
% Arrange the final tensor so that the columns are normalized.
P = arrange(P);

if printitn>0
  normresidual = sqrt( normX^2 + norm(P)^2 - 2 * innerprod(X,P) );
  fit = 1 - (normresidual / normX); %fraction explained by model
  fprintf(' Final fit = %e \n', fit);
end

return;

%% 
function x = setparam(opts,name,default)
if isfield(opts,name);
    x = opts.(name);
else
    x = default;
end
