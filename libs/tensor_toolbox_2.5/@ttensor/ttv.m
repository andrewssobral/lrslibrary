function c = ttv(a,v,dims)
%TTV Tensor times vector for ttensor.
%
%   Y = TTV(X,A,N) computes the product of ttensor X with a
%   (column) vector A.  The integer N specifies the dimension in X
%   along which A is multiplied.  If size(A) = [I,1], then X must have
%   size(X,N) = I.  Note that ndims(Y) = ndims(X) - 1 because the N-th
%   dimension is removed.
%
%   Y = TTV(X,{A1,A2,...}) computes the product of ttensor X with a
%   sequence of vectors in the cell array.  The products are computed
%   sequentially along all dimensions (or modes) of X. The cell array
%   contains ndims(X) vectors.
%
%   Y = TTV(X,{A1,A2,...},DIMS) computes the sequence tensor-vector
%   products along the dimensions specified by DIMS.
%
%   See also TENSOR/TTV, TTENSOR, TTENSOR/TTM.
%
%MATLAB Tensor Toolbox.
%Copyright 2010, Sandia Corporation. 

% This is the MATLAB Tensor Toolbox by Brett Bader and Tamara Kolda. 
% http://csmr.ca.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2010) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in tensor_toolbox/LICENSE.txt
% $Id: ttv.m,v 1.6 2010/03/19 23:46:31 tgkolda Exp $


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

% Create w to be multiplied with a.core
w = cell(ndims(a),1);
for i = 1:numel(dims)
    w{dims(i)} = a.u{dims(i)}' * v{vidx(i)};
end

% Create new core
newcore = ttv(a.core,w,dims);

% Create final result
if isempty(remdims)
    c = newcore;
else
    c = ttensor(newcore,a.u{remdims});
end
