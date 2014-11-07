function Y = ttm_me(X, U, edims, sdims, tflag)
%TTM_ME Memory-efficient sptensor times matrix.
%
%   Y = TTM_ME(X, U, EDIMS, SDIMS, TFLAG) handles some dimensions
%   elementwise and others in the standard way. Here, X is a sparse tensor
%   (sptensor), U is a cell array of matrices of length ndims(X), 
%   EDIMS specifies the dimensions that are to be handled elementwise, and
%   SDIMS specifies those dimensions that are to be handled in the standard
%   way, and TFLAG indicates to multiply with matrix transpose. Ultimately,
%   the result should be equivalent to called TTM(X, U, [EDIMS SDIMS]) but
%   will use less memory in the computation.
%
%   See also TENSOR_TOOLBOX, SPTENSOR/TTM, TUCKER_ME, TTM_ME_PARTITION.
%
%   Code by Tamara Kolda and Jimeng Sun, 2008. Adapted with permission from
%   the MATLAB Tensor Toolbox function @sptensor/ttm.
%
%   Based on the paper:
%   T. G. Kolda and J. Sun. Scalable Tensor Decompositions for Multi-aspect
%   Data Mining. In: ICDM 2008: Proceedings of the 8th IEEE International
%   Conference on Data Mining, December 2008.  

% $Id: ttm_me.m,v 1.2 2009/07/08 00:03:52 tgkolda Exp $

%% Setup and error checking

% Check number of input arguments
if (nargin < 4)
    error('TTM_ME requires four arguments');
end
if (nargin==4)
    tflag = '';
else
    tflag = 't';
end
% Check that X is a sparse tensor
if ~isa(X, 'sptensor')
    error('Input tensor X must be sparse');
end

% Set number of dimensions of X
N = ndims(X);

% Check that U is a cell array
if ~iscell(U)
    error('U must be a cell array');
end

% Check that the cell array U
if numel(U) ~= N
    error('Incorrect number of elements in U');
end

% Check that every member of edims is in 1:N
tf = ismember(edims, 1:N);
if min(tf) == 0
    error('Invalid dimensions specified');
end

% Check that every member of sdims is in 1:N
tf = ismember(sdims, 1:N);
if min(tf) == 0
    error('Invalid dimensions specified');
end

% Check that edims and sdims are distinct
idx = intersect(edims, sdims);
if ~isempty(idx)
    error('Invalid dimensions specified');
end

%% Check for special case of empty edims
% This means that none of the modes is to be handled elementwise
if isempty(edims)
    Y = ttm(X, U, sdims, tflag);
    return;
end

%% Determine some sizes and set up some variables and mappings

% Determine size of Y (final result)
sizY = size(X);
for n = union(edims, sdims)
    if strcmp(tflag,'t')
        sizY(n) = size(U{n}, 2);
    else
        sizY(n) = size(U{n}, 1);
    end
end

% Allocate space for Y (final result)
Y = tenzeros(sizY);

% Set up cell array of vectors for elementwise computations
v = cell(length(edims), 1);

% Set up mapping from sdims on X to appropriate dimensions on Z
% (zdims) as well as array of matrices (M).
zmap = -1 * ones(1,N);
j = 0;
for i = 1:N
    if ~ismember(i, edims)
        j = j+1;
        zmap(i) = j;
        M{j} = U{i};
    end
end
zdims = zmap(sdims);

%% Main Loop
for i = 1:prod(sizY(edims))
    % Get the subscripts of the rows to be extracted
    rsubs = tt_ind2sub(sizY(edims), i);  

    % Extract the appropriate rows of the U matrices
    for j = 1:length(edims)
        if strcmp(tflag,'t')
            v{j} = U{edims(j)}(:, rsubs(j));  
        else
            v{j} = (U{edims(j)}(rsubs(j), :))';
        end
    end
    
    % Create argument to pass into tensor/subsasgn
    rsubsarg = makesubsarg(N, edims, rsubs);
         
    % Assign to appropriate part of Y
    if isempty(zdims)
        % Case 1: Assigning a single element of Y
        t1 = ttv(X, v, edims);
        Y = subsasgn(Y, rsubsarg, t1);
    else
        % Case 2: Assigning an entire subtensor
        Z = ttv(X, v, edims);
        if nnz(Z)==0
            tmp = 0;
        else
            tmp = full(ttm(Z, M, zdims,tflag));
        end            
        Y = subsasgn(Y, rsubsarg, tmp);
    end
    
end

%%
function s = makesubsarg(n, edims, rsubs)
s.type = '()';
s.subs = cell(1, n);
j = 1;
for i = 1:n
    if ismember(i, edims)
        s.subs{i} = [rsubs(j)];
        j = j + 1;
    else
        s.subs{i} = ':';
    end
end