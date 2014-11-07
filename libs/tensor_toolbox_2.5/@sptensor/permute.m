function t = permute(t,order)
%PERMUTE Rearrange the dimensions of a sparse tensor.
%
%   B = PERMUTE(A,ORDER) rearranges the dimensions of A so that they
%   are in the order specified by the vector ORDER. The result has the
%   same values of A, but the order of the subscripts needed to access
%   any particular element are rearranged as specified by ORDER.
%
%   See also SPTENSOR, PERMUTE.
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


% Error checking
if (ndims(order) ~= 2) || (size(order,1) ~= 1) 
    error('ORDER must be a row vector');
end
   
% Check that the permuation is valid
if ~isequal(sort(order),1:ndims(t))
    error('Invalid permutation.');
end

% Do the permutation
if ~isempty(t.subs)
    t.subs = t.subs(:,order);
end
t.size = t.size(order);
