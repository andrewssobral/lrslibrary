function X = ttm(X,V,varargin)
%TTM Tensor times matrix for ttensor.
%
%   Y = TTM(X,A,N) computes the n-mode product of the ttensor X with a
%   matrix A; i.e., X x_N A.  The integer N specifies the dimension
%   (or mode) of X along which A should be multiplied.  If size(A) =
%   [J,I], then X must have size(X,N) = I.  The result will be a
%   ttensor of the same order and size as X except that size(Y,N) = J.
%
%   Y = TTM(X,{A,B,C,...}) computes the n-mode product of the ttensor
%   X with a sequence of matrices in the cell array.  The n-mode
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
%   X = ttensor(tensor(rand(2,2,2,2)),{rand(5,2),rand(3,2),rand(4,2),rand(2,2)});
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
%   See also TTENSOR, TTENSOR/ARRANGE, TENSOR/TTM
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
% $Id: ttm.m,v 1.7 2010/03/19 23:46:31 tgkolda Exp $


%%%%%%%%%%%%%%%%%%%%%%
%%% ERROR CHECKING %%%
%%%%%%%%%%%%%%%%%%%%%%

% Check the number of arguments
if (nargin < 2)
    error('TTM requires at least two arguments.');
end

% Check for transpose option
isTranspose = false;
if numel(varargin) > 0
  if isnumeric(varargin{1});
    dims = varargin{1};
  end
  isTranspose =  (ischar(varargin{end}) && (varargin{end} == 't'));
end

% Check for dims argument
if ~exist('dims','var')
    dims = [];
end

% Check that 2nd argument is cell array. If not, recall with V as a
% cell array with one element.
if ~iscell(V)
    X = ttm(X,{V},dims,varargin{end});
    return;
end

% Get sorted dims and index for multiplicands
[dims,vidx] = tt_dimscheck(dims,ndims(X),numel(V));

% Determine correct size index
if isTranspose
  j = 1; 
else
  j = 2;
end

% Check that each multiplicand is the right size.
for i = 1:numel(dims)
    if (ndims(V) ~= 2) || (size(V{vidx(i)},j) ~= size(X,dims(i)))
        error('Multiplicand is wrong size');
    end
end

% Do the multiplications in the specified modes. 
for i = 1:numel(dims) 
  if isTranspose
    X.u{dims(i)} = V{vidx(i)}'* X.u{dims(i)};
  else
    X.u{dims(i)} = V{vidx(i)} * X.u{dims(i)};
  end
end
