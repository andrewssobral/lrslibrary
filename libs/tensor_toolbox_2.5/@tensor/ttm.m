function Y = ttm(X,V,varargin)
%TTM Tensor times matrix.
%
%   Y = TTM(X,A,N) computes the n-mode product of tensor X with a
%   matrix A; i.e., X x_N A.  The integer N specifies the dimension
%   (or mode) of X along which A should be multiplied.  If size(A) =
%   [J,I], then X must have size(X,N) = I.  The result will be the
%   same order and size as X except that size(Y,N) = J.
%
%   Y = TTM(X,{A,B,C,...}) computes the n-mode product of tensor X
%   with a sequence of matrices in the cell array.  The n-mode
%   products are computed sequentially along all dimensions (or modes)
%   of X. The cell array contains ndims(X) matrices.
%
%   Y = TTM(X,{A,B,C,...},DIMS) computes the sequence tensor-matrix
%   products along the dimensions specified by DIMS.
%
%   Y = TTM(...,'t') performs the same computations as above except
%   the matrices are transposed.
%
%   Examples
%   X = tensor(rand(5,3,4,2));
%   A = rand(4,5); B = rand(4,3); C = rand(3,4); D = rand(3,2);
%   Y = ttm(X, A, 1)         %<-- computes X times A in mode-1
%   Y = ttm(X, {A,B,C,D}, 1) %<-- same as above
%   Y = ttm(X, A', 1, 't')   %<-- same as above
%   Y = ttm(X, {A,B,C,D}, [1 2 3 4]) %<-- 4-way multiply
%   Y = ttm(X, {D,C,B,A}, [4 3 2 1]) %<-- same as above
%   Y = ttm(X, {A,B,C,D})            %<-- same as above
%   Y = ttm(X, {A',B',C',D'}, 't')   %<-- same as above
%   Y = ttm(X, {C,D}, [3 4])     %<-- X times C in mode-3 & D in mode-4
%   Y = ttm(X, {A,B,C,D}, [3 4]) %<-- same as above
%   Y = ttm(X, {A,B,D}, [1 2 4])   %<-- 3-way multiply
%   Y = ttm(X, {A,B,C,D}, [1 2 4]) %<-- same as above
%   Y = ttm(X, {A,B,D}, -3)        %<-- same as above
%   Y = ttm(X, {A,B,C,D}, -3)      %<-- same as above
%
%   See also TENSOR, TENSOR/TTT, TENSOR/TTV.
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



%% Check the number of arguments
if (nargin < 2)
    error('TTM requires at least two arguments.');
end

%% Create 'n' and 'tflag' arguments from varargin
n = 1:ndims(X);
tflag = '';
if numel(varargin) == 1
    if ischar(varargin{1})
        tflag = varargin{1};
    else
        n = varargin{1};
    end
elseif numel(varargin) == 2
    n = varargin{1};
    tflag = varargin{2};
end

%% Handle cell array
if iscell(V)   
    % Copy n into dims
    dims = n;
    % Check that the dimensions are valid
    [dims,vidx] = tt_dimscheck(dims,ndims(X),numel(V));
    % Calculate individual products
    Y = ttm(X, V{vidx(1)}, dims(1), tflag);
    for i = 2 : numel(dims)
        Y = ttm(Y, V{vidx(i)}, dims(i), tflag);
    end
    % All done
    return;
end

%% Check the second argument
if ndims(V) ~= 2
    error('tensor/ttm: 2nd argument must be a matrix.');
end

%% Check n
if (numel(n) ~= 1 || (n < 0) || n > ndims(X))
    error('Dimension N must be between 1 and NDIMS(X).');
end

%% COMPUTE THE PRODUCT 

N = ndims(X);
sz = size(X);
order = [n,1:n-1,n+1:N];
newdata = double(permute(X,order));
newdata = reshape(newdata,sz(n),prod(sz([1:n-1,n+1:N])));
if tflag == 't'
    newdata = V'*newdata;
    p = size(V,2);
else
    newdata = V*newdata;
    p = size(V,1);
end
newsz = [p,sz(1:n-1),sz(n+1:N)];
Y = tensor(newdata,newsz);
Y = ipermute(Y,order);


return;

end
