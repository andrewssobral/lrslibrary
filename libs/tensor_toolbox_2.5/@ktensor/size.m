function m = size(t,idx)
%SIZE Size of ktensor.
%  
%   D = SIZE(T) returns the size of the tensor. 
%
%   I = SIZE(T,DIM) returns the size of the dimension specified by
%   the scalar DIM.
%
%   See also KTENSOR, KTENSOR/NDIMS.
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


if isempty(t.lambda)
    m = [];
end

if exist('idx','var')
    m = size(t.u{idx}, 1);
else
    for i = 1 : ndims(t)
	m(i) = size(t.u{i}, 1);
    end
end
