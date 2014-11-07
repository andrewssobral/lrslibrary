function a = sptenmat(varargin)
%SPTENMAT Matricized sparse tensor stored as a sparse 2D array.
%
%   A = SPTENMAT(T, RDIMS) creates a sparse matrix representation of
%   an sptensor T.  The dimensions (or modes) specified in RDIMS map
%   to the rows of the matrix, and the remaining dimensions (in
%   ascending order) map to the columns.
%
%   A = SPTENMAT(T, CDIMS, 't') does the same as above, but instead
%   the column dimensions are specified, and the remaining dimensions
%   (in ascending order) map to the rows.
%
%   A = SPTENMAT(T, RDIMS, CDIMS) creates a sparse matrix
%   representation of sptensor T.  The dimensions specified in RDIMS
%   map to the rows of the matrix, and the dimensions specified in
%   CDIMS map to the columns, in the order given.
%
%   A = SPTENMAT(T, RDIM, STR) creates the same matrix representation
%   as above, except only one dimension in RDIM maps to the rows of
%   the matrix, and the remaining dimensions span the columns in an
%   order specified by the string argument STR as follows:
%
%     'fc' - Forward cyclic.  Order the remaining dimensions in the
%            columns by [RDIM+1:ndims(T), 1:RDIM-1].  This is the
%            ordering defined by Kiers.
%
%     'bc' - Backward cyclic.  Order the remaining dimensions in the
%            columns by [RDIM-1:-1:1, ndims(T):-1:RDIM+1].  This is the
%            ordering defined by De Lathauwer, De Moor, and Vandewalle.
%
%   A = SPTENAMT(B,RDIMS,CDIMS,TSIZE) creates a sptenmat from a matrix B
%   along with the mappings of the row (RDIMS) and column indices (CDIMS)
%   and the size of the original tensor (TSIZE). 
%
%   A = SPTENMAT(SUBS, VALS, RDIMS, CDIMS, TSIZE) creates a sptenmat
%   from a set of 2D subscripts (SUBS) and values (VALS) along with
%   the mappings of the row (RDIMS) and column indices (CDIMS) and the
%   size of the original tensor (TSIZE).
%
%   A = SPTENMAT is the empty constructor.
%
%   See also SPTENSOR, TENMAT.
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


%----------
% EMPTY CONSTRUCTOR
%----------
if (nargin == 0)
    a.tsize = [];
    a.rdims = [];
    a.cdims = [];
    a.subs = [];
    a.vals = [];
    a = class(a, 'sptenmat');
    return;
end

%----------
% COPY CONSTRUCTOR
%----------
if (nargin == 1) && isa(varargin{1},'sptenmat')
    t = varargin{1};
    a.tsize = t.tsize;
    a.rdims = t.rdims;
    a.cdims = t.cdims;
    a.subs = t.subs;
    a.vals = t.vals;
    a = class(a, 'sptenmat');
    return;
end

%----------
% CONVERT LIST OF SUBS/VALS 
%----------
if (nargin == 5)

    subs = varargin{1};
    vals = varargin{2};
    rdims = varargin{3};
    cdims = varargin{4};
    tsize = varargin{5};
    
    % Error check
    n = numel(tsize);
    if ~isequal(1:n, sort([rdims cdims]))
        error('Incorrect specification of dimensions');
    elseif ~isempty(subs) && prod(tsize(rdims)) < max(subs(:,1))
        error('Invalid row index');
    elseif ~isempty(subs) && prod(tsize(cdims)) < max(subs(:,2))
        error('Invalid column index');
    end

    % Sum any duplicates
    if isempty(subs)
        newsubs = [];
        newvals = [];
    else
        % Identify only the unique indices
        [newsubs,junk,loc] = unique(subs,'rows');

        % Sum the corresponding values
        newvals = accumarray(loc,vals,[size(newsubs,1) 1]);
    end
    
    % Find the nonzero indices of the new values
    nzidx = find(newvals);
    newsubs = newsubs(nzidx,:);
    newvals = newvals(nzidx);
  
    % Save class variables
    a.tsize = tsize;
    a.rdims = rdims;
    a.cdims = cdims;
    a.subs = newsubs;
    a.vals = newvals;
    a = class(a, 'sptenmat');
    return;

end

%----------
% CONVERT SPARSE or DENSE MATLAB MATRIX
%----------
if (nargin == 4)

    B = varargin{1};
    [i,j,vals] = find(B);
    subs = [i j];
    rdims = varargin{2};
    cdims = varargin{3};
    tsize = varargin{4};
    
    % Error check
    n = numel(tsize);
    if ~isequal(1:n, sort([rdims cdims]))
        error('Incorrect specification of dimensions');
    elseif ~isempty(subs) && prod(tsize(rdims)) < max(subs(:,1))
        error('Invalid row index');
    elseif ~isempty(subs) && prod(tsize(cdims)) < max(subs(:,2))
        error('Invalid column index');
    end

    % Save class variables
    a.tsize = tsize;
    a.rdims = rdims;
    a.cdims = cdims;
    a.subs = subs;
    a.vals = vals;
    a = class(a, 'sptenmat');
    return;

end


%----------
% CONVERT SPTENSOR 
%----------

if (nargin < 2)  ||  (nargin > 3)
  error('Incorrect number of arguments.');
end

% Save the size of T and the number of dimensions
T = varargin{1};
tsize = size(T);
tsubs = T.subs;
tvals = T.vals;
n = ndims(T);

% Figure out which dimensions get mapped where
if (nargin == 2) 
    rdims = varargin{2};
    cdims = setdiff(1:n, rdims);
elseif isa(varargin{3},'char')
    switch varargin{3}
        case 't'                        % Transpose
            cdims = varargin{2};
            rdims = setdiff(1:n, cdims);
        case 'fc'                       % Forward cyclic
            rdims = varargin{2};
            if (numel(rdims) ~= 1)
                error('Only one row dimension if third argument is ''fc''.');
            end
            cdims = [rdims+1:n, 1:rdims-1];
        case 'bc'                       % Backward cyclic
            rdims = varargin{2};
            if (numel(rdims) ~= 1)
                error('Only one row dimension if third argument is ''bc''.');
            end
            cdims = [rdims-1:-1:1, n:-1:rdims+1];
        otherwise
            error('Unrecognized option');
    end
else
    rdims = varargin{2};
    cdims = varargin{3};
end

% Error check
if ~isequal(1:n, sort([rdims cdims]))
    error('Incorrect specification of dimensions');
end

% Extract the appropriate sizes
rsize = tsize(rdims);
csize = tsize(cdims);

% Reshape by transforming the indices
if isempty(rsize)
    ridx = ones(nnz(T),1);
elseif isempty(tsubs)
    ridx = [];
else
    ridx = tt_sub2ind(rsize,tsubs(:,rdims));
end

if isempty(csize)
    cidx = ones(nnz(T),1);
elseif isempty(tsubs)
    cidx = [];
else
    cidx = tt_sub2ind(csize,tsubs(:,cdims));
end

% Save class variables
a.tsize = tsize;
a.rdims = rdims;
a.cdims = cdims;
a.subs = [ridx, cidx];
a.vals = tvals;
a = class(a, 'sptenmat');

    
