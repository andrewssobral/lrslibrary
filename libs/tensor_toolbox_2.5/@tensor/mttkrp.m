function V = mttkrp(X,U,n)
%MTTKRP Matricized tensor times Khatri-Rao product for tensor.
%
%   V = MTTKRP(X,U,n) efficiently calculates the matrix product of the
%   n-mode matricization of X with the Khatri-Rao product of all
%   entries in U, a cell array of matrices, except the nth.  How to
%   most efficiently do this computation depends on the type of tensor
%   involved.
%
%   See also TENSOR, TENMAT, KHATRIRAO
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


N = ndims(X);
if (N < 2)
    error('MTTKRP is invalid for tensors with fewer than 2 dimensions');
end

if (length(U) ~= N)
    error('Cell array is the wrong length');
end

if n == 1
    R = size(U{2},2);
else
    R = size(U{1},2);
end

for i = 1:N
   if i == n, continue; end
   if (size(U{i},1) ~= size(X,i)) || (size(U{i},2) ~= R)
       error('Entry %d of cell array is wrong size', i);
   end
end

Xn = permute(X,[n 1:n-1,n+1:N]);
Xn = reshape(Xn.data, size(X,n), prod(size(X))/size(X,n)); %#ok<PSIZE>
Z = khatrirao(U{[1:n-1,n+1:N]},'r');
V = Xn*Z;
