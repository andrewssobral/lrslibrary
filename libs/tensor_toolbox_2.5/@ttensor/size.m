function m = size(t,idx)
%SIZE Size of a ttensor.
%  
%   D = SIZE(T) returns the size of the tensor. 
%
%   I = size(T,DIM) returns the size of the dimension specified by
%   the scalar DIM.
%
%   See also TTENSOR, TTENSOR/NDIMS.
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
% $Id: size.m,v 1.8 2010/03/19 23:46:31 tgkolda Exp $

if ndims(t) == 0
    m = [];
end

if exist('idx','var')
    m = size(t.u{idx}, 1);
else
    for i = 1 : ndims(t)
	m(i) = size(t.u{i}, 1);
    end
end
