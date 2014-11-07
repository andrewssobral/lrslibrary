function m = size(t,idx)
%SIZE Tensor dimensions.
%  
%   D = SIZE(T) returns the sizes of each dimension of tensor X in a
%   vector D with ndims(X) elements.
%
%   I = size(T,DIM) returns the size of the dimension specified by
%   the scalar DIM.
%
%   Examples
%   A = rand(3,4,2,1); T = tensor(A,[3 4 2 1]);
%   size(A) %<-- returns a length-3 vector
%   size(T) %<-- returns a length-4 vector
%   size(A,2) %<-- returns 4
%   size(T,2) %<-- same
%   size(A,5) %<-- returns 1
%   size(T,5) %<-- ERROR!
%
%   See also TENSOR, TENSOR/NDIMS, SIZE.
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



if exist('idx','var')
    m = t.size(idx);
else
    m = t.size;
end
