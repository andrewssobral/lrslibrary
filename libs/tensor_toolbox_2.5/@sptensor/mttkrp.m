function V = mttkrp(X,U,n)
%MTTKRP Matricized tensor times Khatri-Rao product for sparse tensor.
%
%   V = MTTKRP(X,U,n) efficiently calculates the matrix product of the
%   n-mode matricization of X with the Khatri-Rao product of all
%   entries in U, a cell array of matrices, except the nth.  How to
%   most efficiently do this computation depends on the type of tensor
%   involved.
%
%   See also SPTENSOR, TENSOR/MTTKRP, SPTENSOR/TTV
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


% In the sparse case, it is most efficient to do a series of TTV operations
% rather than forming the Khatri-Rao product.

N = ndims(X);

if (n == 1)
    R = size(U{2},2);
else
    R = size(U{1},2);
end

V = zeros(size(X,n),R);
for r = 1:R
    % Set up cell array with appropriate vectors for ttv multiplication
    Z = cell(N,1);
    for i = [1:n-1,n+1:N]
        Z{i} = U{i}(:,r);
    end
    % Perform ttv multiplication
    V(:,r) = double(ttv(X, Z, -n));
end
