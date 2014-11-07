function e = end(X,k,n)
%END Last index of indexing expression for ktensor.
%
%   The expression X(end,:,:) will call END(X,1,3) to determine
%   the value of the first index.
%
%   See also KTENSOR, KTENSOR/SUBSREF, END.
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


%TODO (after 2.0 release): Resolve ambiguity w.r.t X{end}and X(end,1,1)
%for 1st-order tensors.

if n > ndims(X)
  error('Subscript out of range.');
end

if (n ~= 1)
  e = size(X,k);
else
  e = ndims(X);
end
