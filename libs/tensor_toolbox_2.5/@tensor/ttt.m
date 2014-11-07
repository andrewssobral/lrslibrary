function c = ttt(varargin)
%TTT Tensor mulitplication (tensor times tensor).
% 
%   TTT(X,Y) computes the outer product of tensors X and Y.
%
%   TTT(X,Y,XDIMS,YDIMS) computes the contracted product of tensors 
%   X and Y in the dimensions specified by the row vectors XDIMS and 
%   YDIMS.  The sizes of the dimensions specified by XDIMS and YDIMS 
%   must match; that is, size(X,XDIMS) must equal size(Y,YDIMS). 
%
%   TTT(X,Y,DIMS) computes the inner product of tensors X and Y in the
%   dimensions specified by the vector DIMS.  The sizes of the
%   dimensions specified by DIMS must match; that is, size(X,DIMS) must
%   equal size(Y,DIMS). 
%
%   Examples
%   X = tensor(rand(4,2,3));
%   Y = tensor(rand(3,4,2));
%   Z = ttt(X,Y) %<-- outer product of X and Y
%   Z = ttt(X,X,1:3) %<-- inner product of X with itself
%   Z = ttt(X,Y,[1 2 3],[2 3 1]) %<-- inner product of X & Y
%   Z = ttt(X,Y,[1 3],[2 1]) %<-- product of X & Y along specified dims
%
%   See also TENSOR, TENSOR/TTM, TENSOR/TTV.
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
    error('TTT requires at least two arguments.');
end

% Check the first argument
if ~isa(varargin{1}, 'tensor')
    error('First argument must be a tensor.');
else
    a = varargin{1};
end

% Check the second argument
if ~isa(varargin{2}, 'tensor')
    error('Second argument must be a tensor.');
else
    b = varargin{2};
end

% Optional 3rd argument
if nargin >= 3
    adims = varargin{3};
else
    adims = [];
end

% Optional 4th argument
if nargin >= 4
    bdims = varargin{4};
else
    bdims = adims;
end

if ~isequal(size(a,adims),size(b,bdims))
    error('Specified dimensions do not match.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPUTE THE PRODUCT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Avoid transpose by reshaping A and computing C = A * B
amatrix = tenmat(a,adims,'t');
bmatrix = tenmat(b,bdims);
cmatrix = amatrix * bmatrix;

% Check whether or not the result is a scalar.
if isa(cmatrix,'tenmat')
    c = tensor(cmatrix);
else
    c = cmatrix;
end
