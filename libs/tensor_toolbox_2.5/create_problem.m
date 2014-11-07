function [info,params] = create_problem(varargin)
%CREATE_PROBLEM Create test problems for tensor factorizations.
%
%   INFO = CREATE_PROBLEM('Param',value,...) creates a tensor factorization
%   test problem. It generates a solution corresponding to a ktensor or a
%   ttensor, and then it generates an example data tensor which has that
%   underlying factorization. The data tensor can be optionally corrupted
%   with noise, generated via a special procedure to produce a sparse
%   tensor, or even have missing data.
%
%   [INFO,PARAMS] = CREATE_PROBLEM(...) also returns the parameters that
%   were used to generate the problem. An indentical problem can be
%   generated via INFO_COPY = CREATE_PROBLEM(PARAMS).
%  
%   --- General ---
%
%   'State' - State of the random number generator. This can be used
%   to reproduce results.
%
%   --- Solution Parameters ---
%
%   The desired solution will be returned in the "Soln" field of INFO. It
%   will be either a ktensor or ttensor, depending on the user option
%   'Type', described below. 
%
%   'Soln' - A CP or Tucker tensor that is the desired solution. Renders
%   all other solution parameters obsolete. Default: []
%
%   'Type' - CP or Tucker. Default: 'CP'
%
%   'Size' - Size of the tensor. Default: [100 100 100]
%
%   'Num_Factors' - Number of factors. Default: 2 (or [2 2 2] for Tucker)
%
%   'Symmetric' - List of modes that should be symmetric, i.e., have
%   identical factor matrices. This cell array of lists of modes if
%   different subsets should be symmetric, i.e., {[1 2],[3 4]} says that
%   modes 1 and 2 and 3 and 4 are symmetric.
%
%   'Factor_Generator' - Method to be used to generate the factor matrices.
%   Options are 'rand' (uniform on [0,1]), 'randn' (standard normal
%   distribution), 'orthogonal', or 'stochastic' (uniform on [0,1]
%   with column sums rescaled to 1). Alternatively, pass in a function that
%   accepts two arguments (the size of the matrix) and generates the
%   desired factor. Default: 'randn'
%
%   'Lambda_Generator' - Method used to genate the lambda vector for CP
%   solutions. The choices are the same as for the 'Factor_Generator'.
%   Default: 'rand'
%
%   'Core_Generator' - Method used to generate the core tensor for Tucker
%   solutions. The choices are 'rand' and 'randn' (as described above).
%   Alternatively, pass in a function that accepts a vector-valued size and
%   generates a tensor of the specified size. Default: 'randn'
%
%   --- Missing Data Parameters ---
%
%   The missing data pattern will be return in the "Pattern" data
%   field of INFO. If there is no missing data, then this will just be an
%   empty array. Otherwise, it will be a tensor that is zero whereever data
%   is missing and one elsewhere. 
%
%   'M' - The proportion of missing data *or* a tensor or sptensor that
%   contains the missing data pattern as described above. Default: 0 
%
%   'Sparse_M' - Generate sparse rather than dense missing data pattern
%   tensor. Only useful for large tensors that don't easily fit in memory
%   and when M > 80%. Default: false. 
%
%   --- Data Parameters ---
%
%   The data to be factorized will be returned in the "Data" field of INFO.
%   It will have zeros for any entries that are missing (though not all
%   zeros necessarily correspond to missing data).
%
%   'Sparse_Generation' - Generate a sparse tensor via a special procedure
%   that works only for ktensor's (CP) that can be scaled so that the
%   column factors and lambda are stochastic. Note that this geneartion
%   procedure will modify lambda vector in the solution so that it is
%   appropriately scaled for the number of inserted nonzeros. A value of
%   zero means no sparse generation, and any positive value is the number
%   of nonzeros to be inserted. Any value in the range (0,1) will be
%   interpreted as a percentage. The procedure is incompatible with missing
%   data. Default: 0 (no sparse geneartion).
%
%   'Noise' - Amount of Gaussian noise to add. Let N be a "noise"
%   tensor with entries drawn from the standard norm distribution, and
%   let Y be the noise-free tensor, i.e. Y = full(K). Then Z = Y + eta
%   * norm(Y,'fro') / norm(N,'fro') * N is the noisy version of the
%   tensor where eta is the percentage of noise to add. If the data tensor
%   is sparse (either due to sparse generation or sparsity due to missing
%   data), then noise is only generated at the nonzero entries. 
%   Default: 0.10
%
%   Examples:
%   % Create a 100 x 100 x 100 problem with 5 factors (each entry from the
%   % standard normal distribution) and 10% noise with diagonal lambda
%   % values of all ones.
%   info = createcreate_problem('Lambda_Generator', @ones);
%
%   % Same as above except that the we use a special function to generate
%   % factor matrices with a constant congruence of 0.9.
%   info = create_problem('Factor_Generator', @(m,n) tt_ccong(m,n,.9), ...
%     'Lambda_Generator', @ones);
%
%   See also CREATE_GUESS.
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
p.addParamValue('State', defaultStream.State, @(x) true);
p.addParamValue('Soln', [], @(x) isempty(x) || isa(x,'ktensor') || isa(x,'ttensor'));
p.addParamValue('Type', 'CP', @(x) ismember(lower(x),{'cp','tucker'}));
p.addParamValue('Size', [10 10 10], @all);
p.addParamValue('Num_Factors', 2, @all);
p.addParamValue('Factor_Generator', 'randn', @is_valid_matrix_generator);
p.addParamValue('Lambda_Generator', 'rand', @is_valid_matrix_generator);
p.addParamValue('Core_Generator', 'randn', @is_valid_tensor_generator);
p.addParamValue('M', 0, @(x) is_missing_data(x) || (x == 0));
p.addParamValue('Sparse_M', false, @islogical);
p.addParamValue('Sparse_Generation', 0, @(x) x >= 0);
p.addParamValue('Symmetric', []);
p.addParamValue('Noise', 0.10, @(x) x >= 0 & x < 1);
% p.addParamValue('Rtest', 0, @(x) isscalar(x) & x >= 0);
% p.addParamValue('Init_Type', 'random', @(x) ismember(x,{'random','nvecs'}));
p.parse(varargin{:});
params = p.Results;

%% Initialize random number generator with specified state.
defaultStream.State = params.State;

%% Error checking
if is_missing_data(params.M) && (params.Sparse_Generation > 0)
    error('Cannot combine missing data and sparse generation');
end

if strcmpi(params.Type, 'tucker') && (params.Sparse_Generation > 0)
    error('Sparse generation only supported for CP');
end

%% Check for incompatible parameters
if ~isempty(params.Symmetric)
    if is_missing_data(params.M)
        error('Not set up to generate a symmetric problem with missing data');
    end
    if (params.Sparse_Generation ~= 0)
        error('Not set up to generate a sparse symmetric problem');
    end
end        

%% Create return data structure
info = struct;
info.Soln = generate_solution(params);

if is_missing_data(params.M)
    info.Pattern = generate_missing_pattern(size(info.Soln), params);
    info.Data = generate_data(info.Soln, info.Pattern, params);
elseif (params.Sparse_Generation)
    [info.Data, info.Soln] = generate_data_sparse(info.Soln, params);
else
    info.Data = generate_data(info.Soln, [], params);
end

function D = generate_data(S, W, params)
%GENERATE_DATA Generate CP or Tucker data from a given solution.

sz = size(S);
if isempty(W)   
    Rdm = tensor(randn([sz 1 1]), sz);
    Z = full(S);   
    % Use symmetric noise for a symmetric problem
    if ~isempty(params.Symmetric)
        Rdm = symmetrize(Rdm, params.Symmetric);
    end
else
    if isa(W,'sptensor')
        Rdm = sptensor(W.subs,randn(nnz(W),1),W.size);
        Z = W.*S;
    else
        Rdm = W.*tensor(randn([sz 1 1]), sz);
        Z = W.*full(S);
    end
end

D = Z + params.Noise * norm(Z) * Rdm / norm(Rdm);

% Make sure the final result is *absolutely* symmetric
if ~isempty(params.Symmetric)
    D = symmetrize(D, params.Symmetric);
end

function output = prosample(nsamples, prob)
%PROSAMPLE Proportional sampling

% Create bins
bins = min([0 cumsum(prob')],1);
bins(end) = 1;

% Create indices
[~, output] = histc(rand(nsamples,1),bins);


function [Z,S] = generate_data_sparse(S,params)
%GENERATE_DATA_SPARSE Generate sparse CP data from a given solution.

% Error check on S
if any(S.lambda < 0)
    error('All lambda values must be nonnegative');
end
if any(cellfun(@(x) any(x(:)<0), S.U))
    error('All factor matrices must be nonnegative');
end
if ~strcmpi(params.Type,'CP')
    error('Only works for CP');
end
if ~isempty(params.Symmetric)
    warning('Symmetry constraints have been ignored');
end

% Convert S to a probability tensor
P = normalize(S,[],1); 
eta = sum(P.lambda);    
P.lambda = P.lambda / eta;

% Determine how many samples per component
nedges = params.Sparse_Generation;
nd = ndims(P);
nc = size(P.U{1},2);
csample = prosample(nedges, P.lambda);
csums = accumarray(csample,1,[nc 1]);

% Determine subscripts for each randomly sampled entry
sz = size(S);
subs = cell(nc,1);
for c = 1:nc
    nsample = csums(c);
    if nsample == 0
        continue;
    end
    subs{c} = zeros(nsample,nd);
    for d = 1:nd
        PU = P.U{d};
        subs{c}(:,d) = prosample(nsample, PU(:,c));        
    end
end

% Assemble final tensor. Note that duplicates are summed.
allsubs = cell2mat(subs);
Z = sptensor(allsubs,1,sz);

% Rescale S so that it is proportional to the number of edges inserted
S = P;
S.lambda = nedges * S.lambda;

function W = generate_missing_pattern(sz, params)
%GENERATE_MISSING_PATTERN Generate a tensor pattern of missing data.

M = params.M;
S = params.Sparse_M;
if isa(M, 'tensor') || isa(M, 'sptensor')
    W = M;
    return;
end
if M == 0
    W = [];
    return;   
end
if (M < 0.8) && S
    warning('Setting sparse to false because there are less than 80% missing elements');
    S = false;
end
W = tt_create_missing_data_pattern(sz, M, S);

function S = generate_solution(params)
%GENERATE_SOLUTION Generate factor matrices and other data for CP or Tucker

if ~isempty(params.Soln)   
    S = params.Soln;
    return;
end

% Get size of final tensor
sz = params.Size;
nd = length(sz);

% Get size of factors
nfactors = params.Num_Factors;
if numel(nfactors) == 1
    nfactors = nfactors * ones(size(sz));
end
if any(size(nfactors) ~= size(sz))
    error('''Num_Factors'' should either be a single value or the same dimensions as ''Size''.');
end

% Create factor matrices
fgfh = get_generator(params.Factor_Generator);
U = cell(nd,1);
for n = 1:nd
    U{n} = fgfh(sz(n), nfactors(n));
end

if ~isempty(params.Symmetric)
    if ~iscell(params.Symmetric)
        params.Symmetric = {params.Symmetric};
    end
    for i = 1:length(params.Symmetric)
        grp = params.Symmetric{i};
        for j = 2:length(grp)
            U{grp(j)} = U{grp(1)};
        end
    end
end

% Create final ktensor or ttensor
switch lower(params.Type)
    case {'cp'}
        lgfh = get_generator(params.Lambda_Generator);
        lambda = lgfh(nfactors(1),1);
        S = ktensor(lambda,U);
    case {'tucker'}
        cgfh = get_generator(params.Core_Generator);
        core = tensor(cgfh(nfactors));
        S = ttensor(core,U);
    otherwise
        error('Invalid choice for ''Type''');
end

function b = is_valid_matrix_generator(x)
b = isa(x,'function_handle') || ...
    ismember(lower(x),{'rand','randn','orthogonal','stochastic'});

function b = is_valid_tensor_generator(x)
b = isa(x,'function_handle') || ismember(lower(x),{'rand','randn'});

function fh = get_generator(x)
if isa(x,'function_handle')
    fh = x;
    return;
end
switch lower(x)
    case {'randn'}
        fh = @randn;
    case {'rand'}
        fh = @rand;
    case {'orthogonal'}
        fh = @rand_orth_mat;
    case {'stochastic'}
        fh = @rand_column_stochastic;
    otherwise
        error('Invalid choice for generator');
end

function Y = rand_column_stochastic(M,N)
X = rand(M,N);
S = sum(X,1);
Y = X * diag(1./S);

function Y = rand_orth_mat(M,N)
X = tt_RandOrthMat(M);
Y = X(:,1:N);

function tf = is_missing_data(x)
tf = isa(x,'tensor') || isa(x,'sptensor') || (isscalar(x) && (x > 0));
