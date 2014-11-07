function y = contract(x,i,j)
%CONTRACT Contract tensor along two dimensions (array trace).
%
%   Y = CONTRACT(X,I,J) contracts the entries of X along dimensions I
%   and J. Contraction is a generalization of matrix trace. In other
%   words, the trace is performed along the two-dimensional slices
%   defined by dimensions I and J. It is possible to implement tensor
%   multiplication as an outer product followed by a contraction.
%
%   Examples
%   X = tensor(rand(4,3,2)); Y = tensor(rand(3,2,4));
%   Z1 = ttt(X,Y,1,3); %<-- Normal tensor multiplication
%   Z2 = contract(ttt(X,Y),1,6); %<-- Outer product + contract
%   norm(Z1-Z2) %<-- Should be zero
%
%   See also TENSOR, TENSOR/TTT.
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


% Error checking
if x.size(i) ~= x.size(j)
    error('Must contract along equally sized dimensions');
end

% Error checking
if i == j
    error('Must contract along two different dimensions');
end

% Easy case - returns a scalar
if ndims(x) == 2
    y = trace(x.data);
    return;
end

% Remaining dimensions after trace
remdims = setdiff(1:ndims(x),[i j]);

% Size for y
newsize = x.size(remdims);

% Total size of remainder
m =  prod(newsize);

% Number of items to add for trace
n = x.size(i);

% Permute trace dimensions to the end
x = permute(x, [remdims i j]);

% Reshape data to be 3D
data = reshape(x.data, m, n, n);

% Add diagonal entries for each slice
newdata = zeros(m,1);
for i = 1:n
    newdata = newdata + data(:,i,i);
end

% Reshape result
if numel(newsize) > 1
    newdata = reshape(newdata,newsize);
end
y = tensor(newdata,newsize);

