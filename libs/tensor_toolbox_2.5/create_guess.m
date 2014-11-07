function U = create_guess(varargin)
%CREATE_GUESS Creates initial guess for CP or Tucker fitting.
%
%  U = CREATE_GUESS('Param',value,...) creates an initial guess at the
%  factor matrices for a CP or Tucker decomposition. The factors can be
%  generated randomly, random orthogonal, etc. If the tensor is provided,
%  it can be alternatively generated via the HO-SVD. 
%
%   --- Parameters ---
%
%   'Factor_Generator' - Method to be used to generate the factor matrices.
%   Options:
%      - 'rand' (uniform on [0,1])
%      - 'randn' (standard normal distribution)
%      - 'orthogonal'
%      - 'stochastic' (uniform on [0,1] with column sums rescaled to 1)
%      - 'nvecs' (HOSVD solution)
%      - 'pertubation' of the true solution
%   Alternatively, pass in a function that accepts two arguments (the size
%   of the matrix) and generates the desired factor. Default: 'rand' 
%
%   'Size' - Size of the tensor. Required to be specified unless 'Data' or
%   'Soln' is given. Default: [] 
%
%   'Num_Factors' - Number of factors (can be either a single value for CP
%   or a vector for Tucker). Required to be specified unless 'Soln' is
%   given. Default: [] 
%
%   'Data' - The actual tensor to be fit. Required if 'nvecs' is the
%   selected Factor Generator. The 'Size' parameter is ignored if this
%   is specified. Default: []
%
%   'Soln' - The actual solution to the problem. Required if 'pertubation'
%   is the selected Factor Generator. The 'Size' and 'Num_Factors'
%   parameters are ignored if this is specified. Default: []
%
%   'Pertubation' - Size of the pertubation is the 'pertubation' option is
%   selected under 'Factor_Generator'. The pertubation is of the form U+p*N
%   where U is the original factor matrix, N is a noise matrix with entries
%   selected for a standard normal distribution, and p is the pertubation
%   parameter times ||U||/||N||. Default: 0.10
%
%   'Skip' - Specifies mode to skip in initial guess generation (this is
%   useful for ALS). Default: 0 (no skipping)
%
%   'State' - State of the random number generator. This can be used
%   to reproduce results.
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


%% Random set-up
defaultStream = RandStream.getDefaultStream;

%% Parse inputs
p = inputParser;
p.addParamValue('Factor_Generator', 'rand', @(x) isa(x,'function_handle') || ...
    ismember(lower(x),{'rand','randn','orthogonal','stochastic','nvecs','pertubation'}));
p.addParamValue('Size', [], @(x) isempty(x) || all(x));
p.addParamValue('Num_Factors', [], @(x) isempty(x) || all(x));
p.addParamValue('Soln', [], @(x) isempty(x) || isa(x,'ktensor') || isa(x,'ttensor'));
p.addParamValue('Data', [], @(x) isempty(x) || isa(x,'tensor') || isa(x,'sptensor'));
p.addParamValue('Pertubation', 0.10, @(x) x >= 0 & x < 1);
p.addParamValue('Skip', 0);
p.addParamValue('State', defaultStream.State, @(x) true);
p.parse(varargin{:});
params = p.Results;

%% Initialize random number generator with specified state.
defaultStream.State = params.State;

%% Determine problem size
if ~isempty(params.Soln)
    sz = size(params.Soln);
elseif ~isempty(params.Data)
    sz = size(params.Data);
else
    sz = params.Size;
end
if isempty(sz)
    error('Size must be specified');
end
nd = length(sz);
modes = setdiff(1:nd,params.Skip);

%% Determine number of factors
if ~isempty(params.Soln)
    nf = zeros(nd,1);
    for n = 1:nd
        nf(n) = size(params.Soln.U{n},2);
    end
else
    nf = params.Num_Factors;
    if length(nf) == 1
        nf = nf * ones(nd,1);
    end
end

%% Create factor matrices
U = cell(nd,1);
if isa(params.Factor_Generator,'function_handle')
    for n = modes
        U{n} = params.Factor_Generator(sz(n), nf(n));
    end
    return;
end

switch(params.Factor_Generator)
    case 'rand'       
        for n = modes
            U{n} = rand(sz(n), nf(n));
        end
    case 'randn'       
        for n = modes
            U{n} = randn(sz(n), nf(n));
        end
    case 'orthogonal'
        for n = modes
            X = tt_RandOrthMat(sz(n));
            U{n} = X(:,1:nf(n));
        end
    case 'stochastic'
        for n = modes
            X = rand(sz(n), nf(n));
            S = sum(X,1);
            U{n} = X * diag(1./S);
        end
    case 'nvecs'
        if isempty(params.Data)
            error('Data required for nvecs initialization');
        end
        for n = modes
            U{n} = nvecs(params.Data,n,nf(n));
        end
    case 'pertubation'
        if isempty(params.Soln)
            error('Soln required for pertubation initialization');
        end
        for n = modes
            X = params.Soln{n};
            N = rand(size(X));          
            p = params.Pertubation * norm(X,'fro') / norm(N,'fro');
            U{n} = X + p * N;
        end
end
    