function c = ttv(a,v,dims)
%TTV Tensor times vector.
%
%   Y = TTV(X,A,N) computes the product of tensor X with a (column)
%   vector A.  The integer N specifies the dimension in X along which
%   A is multiplied.  If size(A) = [I,1], then X must have size(X,N) =
%   I.  Note that ndims(Y) = ndims(X) - 1 because the N-th dimension
%   is removed.
%
%   Y = TTV(X,{A,B,C,...}) computes the product of tensor X with a
%   sequence of vectors in the cell array.  The products are computed
%   sequentially along all dimensions (or modes) of X. The cell array
%   contains ndims(X) vectors.
%
%   Y = TTV(X,{A,B,C,...},DIMS) computes the sequence of tensor-vector
%   products along the dimensions specified by DIMS.
%
%   Examples
%   X = tensor(rand(5,3,4,2));
%   A = rand(5,1); B = rand(3,1); C = rand(4,1); D = rand(2,1);
%   Y = ttv(X, A, 1) %<-- X times A in mode 1
%   Y = ttv(X, {A,B,C,D}, 1) %<-- same as above
%   Y = ttv(X, {A,B,C,D}, [1 2 3 4]) %<-- All-mode multiply
%   Y = ttv(X, {D,C,B,A}, [4 3 2 1]) %<-- same as above
%   Y = ttv(X, {A,B,C,D}) %<-- same as above
%   Y = ttv(X, {C,D}, [3 4]) %<-- X times C in mode-3 & D in mode-4
%   Y = ttv(X, {A,B,C,D}, [3 4]) %<-- same as above
%   Y = ttv(X, {A,B,D}, [1 2 4]) %<-- 3-way mutplication
%   Y = ttv(X, {A,B,C,D}, [1 2 4]) %<-- same as above
%   Y = ttv(X, {A,B,D}, -3) %<-- same as above
%   Y = ttv(X, {A,B,C,D}, -3) %<-- same as above
%
%   See also TENSOR, TENSOR/TTT, TENSOR/TTM.
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


% Check the number of arguments
if (nargin < 2)
    error('TTV requires at least two arguments.');
end

% Check for 3rd argument
if ~exist('dims','var')
    dims = [];
end

% Check that 2nd argument is cell array. If not, recall with v as a
% cell array with one element.
if ~iscell(v)
    c = ttv(a,{v},dims);
    return;
end

% Get sorted dims and index for multiplicands
[dims,vidx] = tt_dimscheck(dims,ndims(a),numel(v));       

% Check that each multiplicand is the right size.
for i = 1:numel(dims)
    if ~isequal(size(v{vidx(i)}),[size(a,dims(i)) 1])
        error('Multiplicand is wrong size');
    end
end

if exist('tensor/ttv_single','file') == 3
    c = a;
    for i = numel(dims) : -1 : 1
        c = ttv_single(c,v{vidx(i)},dims(i));
    end
    return;
end

% Extract the MDA
c = a.data;

% Permute it so that the dimensions we're working with come last
remdims = setdiff(1:ndims(a),dims);
if (ndims(a) > 1)
    c = permute(c,[remdims dims]);
end

% Do each  multiply in sequence, doing the highest index first,
% which is important for vector multiplies.
n = ndims(a);
sz = a.size([remdims dims]);
for i = numel(dims) : -1 : 1
    c = reshape(c,prod(sz(1:n-1)),sz(n));
    c = c * v{vidx(i)};
    n = n-1;
end

% If needed, convert the final result back to a tensor
if (n > 0)
    c = tensor(c,sz(1:n));
end

