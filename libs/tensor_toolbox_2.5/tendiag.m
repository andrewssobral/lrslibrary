function X = tendiag(v,sz)
%TENDIAG Creates a tensor with v on the diagonal.
%
%   TENDIAG(V) creates a tensor with N dimensions, each of size N, where N
%   is the number of elements of V. The elements of V are placed on the
%   superdiagonal.
%
%   TENDIAG(V,SZ) is the same as above but creates a tensor of size SZ. If
%   SZ is not big enough, the tensor will be enlarged to accommodate the
%   elements of V on the superdiagonal.
%
%   Examples
%   X = tendiag([0.1 0.22 0.333]) %<-- creates a 3x3x3 tensor
%
%   See also TENSOR, SPTENDIAG.
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


% Make sure v is a column vector
v = reshape(v,[numel(v) 1]);

N = numel(v);
if ~exist('sz','var')
    sz = repmat(N,1,N);
end

X = tenzeros(sz);
subs = repmat((1:N)', 1, length(sz));
X(subs) = v;
