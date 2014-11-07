function A = tenmat(varargin)
%TENMAT Create a matricized tensor.
%
%   A = TENMAT(T, RDIMS) creates a matrix representation of a tensor
%   T.  The dimensions (or modes) specified in RDIMS map to the rows
%   of the matrix, and the remaining dimensions (in ascending order)
%   map to the columns.
%
%   A = TENMAT(T, CDIMS, 't') does the same as above, but instead the
%   column dimensions are specified, and the remaining dimensions (in
%   ascending order) map to the rows.
%
%   A = TENMAT(T, RDIMS, CDIMS) creates a matrix representation of
%   tensor T.  The dimensions specified in RDIMS map to the rows of
%   the matrix, and the dimensions specified in CDIMS map to the
%   columns, in the order given.
%
%   A = TENMAT(T, RDIM, STR) creates the same matrix representation as
%   above, except only one dimension in RDIM maps to the rows of the
%   matrix, and the remaining dimensions span the columns in an order
%   specified by the string argument STR as follows:
%
%     'fc' - Forward cyclic.  Order the remaining dimensions in the
%            columns by [RDIM+1:ndims(T), 1:RDIM-1].  This is the
%            ordering defined by Kiers.
%
%     'bc' - Backward cyclic.  Order the remaining dimensions in the
%            columns by [RDIM-1:-1:1, ndims(T):-1:RDIM+1].  This is the
%            ordering defined by De Lathauwer, De Moor, and Vandewalle.
%
%   A = TENMAT(A, RDIMS, CDIMS, TSIZE) creates a tenmat from a matrix
%   A along with the mappings of the row (RDIMS) and column indices
%   (CDIMS) and the size of the original tensor (TSIZE).
%
%   A = TENMAT(B) is the copy constructor for B also a tenmat.
%
%   A = TENMAT is the empty constructor.
%
%   See also TENSOR, SPTENMAT.
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

%
% Includes improvements offered by Marcus Brubaker.


%----------
% Case 0a: Empty Contructor
%----------
if (nargin == 0)
    A.tsize = [];
    A.rindices = [];
    A.cindices = [];
    A.data = [];
    A = class(A, 'tenmat');
    return;
end

%----------
% Case 0b: Copy Contructor
%----------
if (nargin == 1)
    B = varargin{1};
    A.tsize = B.tsize;
    A.rindices = B.rindices;
    A.cindices = B.cindices;
    A.data = B.data;
    A = class(A, 'tenmat');
    return;
end

%----------
% Case I: Called to convert a matrix to a tenmat
%----------
if (nargin == 4)

    data = varargin{1};
    if ~isnumeric(data) || (ndims(data) ~= 2)
        error('A must be a matrix.');
    end
    rdims = varargin{2};
    cdims = varargin{3};
    tsize = varargin{4};

    % Error check
    n = numel(tsize);
    if ~isequal(1:n, sort([rdims cdims]))
        error('Incorrect specification of dimensions');
    elseif (prod(tsize(rdims)) ~= size(data,1))
        error('SIZE(A,1) does not match size specified by RDIMS and SIZE.');
    elseif (prod(tsize(cdims)) ~= size(data,2))
        error('SIZE(A,2) does not match size specified by CDIMS and SIZE.');
    end

    % Save class variables
    A.tsize = tsize;
    A.rindices = rdims;
    A.cindices = cdims;
    A.data = data;
    A = class(A, 'tenmat');
    return;

end

%----------
% Case II: Called to convert an MDA to a tenmat --- recall after
% converting MDA to a tensor.
%----------
if isa(varargin{1},'double')
    A = tenmat(tensor(varargin{1}),varargin{2:nargin});
    return;
end

%----------
% Case III: Convert a tensor to a tenmat
%----------

if (nargin < 2)  ||  (nargin > 3)
    error('Incorrect number of arguments.');
end

% Save the size of T and the number of dimensions
T = varargin{1};
tsize = size(T);
n = ndims(T);

% Figure out which dimensions get mapped where
if (nargin == 2)
    rdims = varargin{2};
    tmp = true(1,n); 
    tmp(rdims) = false; 
    cdims = find(tmp);   % i.e., cdims = setdiff(1:n, rdims);
elseif isa(varargin{3},'char')
    switch varargin{3}
        case 't'                        % Transpose
            cdims = varargin{2};
	    tmp = true(1,n); 
	    tmp(cdims) = false; 
	    rdims = find(tmp);   % i.e., rdims = setdiff(1:n, cdims);
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

% Permute T so that the dimensions specified by RDIMS come first
data = reshape(double(permute(T,[rdims cdims])), prod(tsize(rdims)), prod(tsize(cdims)));

% Save class variables
A.tsize = tsize;
A.rindices = rdims;
A.cindices = cdims;
A.data = data;
A = class(A, 'tenmat');
