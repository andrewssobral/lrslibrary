function [subs,vals] = find(t)
%FIND Find subscripts of nonzero elements in a tensor.
%
%   S = FIND(X) returns the subscripts of the nonzero values in X.
%
%   [S,V] = FIND(X) also returns a column vector of the values.
%
%   Examples: 
%   X = tensor(rand(3,4,2));
%   subs = find(X > 0.5) %<-- find subscripts of values greater than 0.5
%   vals = X(subs) %<-- extract the actual values
%
%   See also TENSOR/SUBSREF, TENSOR/SUBSASGN
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


% Find the *linear* indices of the nonzero elements
idx = find(t.data);

% Convert the linear indices to subscripts
subs = tt_ind2sub(t.size,idx);

% Extract the corresponding values and return as a column vector
if nargout > 1
    if isempty(subs)
        vals = [];
    else
        vals = reshape(t.data(idx), length(idx), 1);
    end
end
