function e = end(X,k,n)
%END Last index of indexing expression for ttensor.
%
%   The expression X(end,:,:) will call END(X,1,3) to determine
%   the value of the first index.
%
%   See also TTENSOR, TTENSOR/SUBSREF, END.
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
% $Id: end.m,v 1.7 2010/03/19 23:46:31 tgkolda Exp $

% Note that this only works with {} because () is not supported by
% subsref.
e = ndims(X);
