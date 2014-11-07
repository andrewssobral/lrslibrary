function Y = scale(X,S,dims)
%SCALE Scale along specified dimensions of tensor.
%
%   Y = SCALE(X,S,DIMS) scales the tensor X along the dimension(s)
%   specified in DIMS using the scaling data in S. If DIMS contains
%   only one dimension, then S can be a column vector. Otherwise, S
%   should be a tensor.
%
%   Examples
%   X = tenones([3,4,5]);
%   S = 10 * [1:5]'; Y = scale(X,S,3)
%   S = tensor(10 * [1:5]',5); Y = scale(X,S,3)
%   S = tensor(1:12,[3 4]); Y = scale(X,S,[1 2])
%   S = tensor(1:12,[3 4]); Y = scale(X,S,-3)
%   S = tensor(1:60,[3 4 5]); Y = scale(X,S,1:3)
%
%   See also TENSOR, TENSOR/COLLAPSE.
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


dims = tt_dimscheck(dims,ndims(X));
remdims = setdiff(1:ndims(X),dims);

% Convert to a matrix so that each column of A can be scaled by a
% vectorized version of S.
A = double(tenmat(X,dims,remdims));

switch(class(S))
    case {'tensor'}
        if ~isequal(size(S), X.size(dims))
            error 'Size mismatch';
        end
        % Vectorize S.
        S = double(tenmat(S,1:ndims(S),[]));
    case {'double'}
        if size(S,1) ~= X.size(dims)
            error 'Size mismatch';
        end
    otherwise
        error('Invalid scaling factor');
end

[m,n] = size(A);

% If the size of S is pretty small, we can convert it to a diagonal matrix
% and multiply by A. Otherwise, we scale A column-by-column.
if (m <= n)
    B = diag(S) * A;
else
    B = zeros(size(A));
    for j = 1:n
        B(:,j) = S .* A(:,j);
    end
end

% Convert the matrix B back into a tensor and return.
Y = tensor(tenmat(B,dims,remdims,X.size));

   


