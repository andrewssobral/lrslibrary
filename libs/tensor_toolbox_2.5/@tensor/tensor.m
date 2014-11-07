function t = tensor(varargin)
%TENSOR Create tensor.
%
%   X = TENSOR(A,SIZ) creates a tensor from the multidimensional
%   array A. The SIZ argument specifies the desired shape of A.
%
%   X = TENSOR(A) creates a tensor from the multidimensional array
%   Z, using SIZ = size(A).
%
%   X = TENSOR(S) copies a tensor S.
%
%   X = TENSOR(A) converts an sptensor, ktensor, ttensor, or tenmat object
%   to a tensor.  
%
%   X = TENSOR creates an empty, dense tensor object.
%
%   Examples
%   X = tensor(rand(3,4,2)) %<-- Tensor of size 3 x 4 x 2
%   Y = tensor(rand(3,1),3) %<-- Tensor of size 3
%   Z = tensor(rand(12,1),[3 4 1]) %<-- Tensor of size 3 x 4 x 1
%
%   See also TENMAT, SPTENSOR, KTENSOR, TTENSOR, TENSOR/NDIMS.
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


% EMPTY/DEFAULT CONSTRUCTOR
if nargin == 0
    t.data = [];
    t.size = [];
    t = class(t, 'tensor');
    return;
end

% CONVERSION/COPY CONSTRUCTORS
% Note that we pass through this if/switch statement if the first argument
% is not any of these cases.
if (nargin == 1)
    v = varargin{1};
    switch class(v)
        case 'tensor',   
            % COPY CONSTRUCTOR
            t.data = v.data;
            t.size = v.size;
            t = class(t, 'tensor');
            return;
        case {'ktensor','ttensor','sptensor'},  
            % CONVERSION
            t = full(v);
            return;
        case 'tenmat', 
            % RESHAPE TENSOR-AS-MATRIX
            % Here we just reverse what was done in the tenmat constructor.
            % First we reshape the data to be an MDA, then we un-permute
            % it using ipermute.
            sz = tsize(v);
            order = [v.rdims,v.cdims];
            data = reshape(v.data, [sz(order) 1 1]);
            if numel(order) >= 2
                t.data = ipermute(data,order);
            else
                t.data = data;
            end              
            t.size = sz;
            t = class(t, 'tensor');
            return;
    end
end

% CONVERT A MULTIDIMENSIONAL ARRAY
if (nargin <= 2)

    % Check first argument
    data = varargin{1};
    if ~isa(data,'numeric') && ~isa(data,'logical')
        error('First argument must be a multidimensional array.')
    end

    % Create or check second argument
    if nargin == 1
        siz = size(data);
    else
        siz = varargin{2};
        if ~isempty(siz) && ndims(siz) ~= 2 && size(siz,1) ~= 1
            error('Second argument must be a row vector.');
        end
    end

    % Make sure the number of elements matches what's been specified
    if isempty(siz)
        if ~isempty(data)
            error('Empty tensor cannot contain any elements');
        end
    elseif prod(siz) ~= numel(data)
        error('Size of data does not match specified size of tensor');
    end
    
    % Make sure the data is indeed the right shape
    if ~isempty(data) && ~isempty(siz)
        data = reshape(data,[siz 1 1]);
    end

    % Create the tensor
    t.data = data;
    t.size = siz;
    t = class(t, 'tensor');
    return;

end


error('Unsupported use of function TENSOR.');


