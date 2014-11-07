function c = ttv(a,v,dims)
%TTV Tensor times vector for ktensor.
%
%   Y = TTV(X,A,N) computes the product of Kruskal tensor X with a
%   (column) vector A.  The integer N specifies the dimension in X
%   along which A is multiplied.  If size(A) = [I,1], then X must have
%   size(X,N) = I.  Note that ndims(Y) = ndims(X) - 1 because the N-th
%   dimension is removed.
%
%   Y = TTV(X,{A1,A2,...}) computes the product of tensor X with a
%   sequence of vectors in the cell array.  The products are computed
%   sequentially along all dimensions (or modes) of X. The cell array
%   contains ndims(X) vectors.
%
%   Y = TTV(X,{A1,A2,...},DIMS) computes the sequence tensor-vector
%   products along the dimensions specified by DIMS.
%
%   See also TENSOR/TTV, KTENSOR, KTENSOR/TTM.
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


%%%%%%%%%%%%%%%%%%%%%%
%%% ERROR CHECKING %%%
%%%%%%%%%%%%%%%%%%%%%%

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

% Figure out which dimensions will be left when we're done
remdims = setdiff(1:ndims(a),dims);

% Collapse dimensions that are being multiplied out
newlambda = a.lambda;
for i = 1:numel(dims) 
    newlambda = newlambda .* ( a.u{dims(i)}' * v{vidx(i)} );
end

% Create final result
if isempty(remdims)
    c = sum(newlambda);
else
    c = ktensor(newlambda,a.u{remdims});
end
