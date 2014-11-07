function y = contract(x,i,j)
%CONTRACT Contract sparse tensor along two dimensions (array trace).
%
%   Y = CONTRACT(X,I,J) contracts the entries of X along dimensions I
%   and J. Contraction is a generalization of matrix trace. In other
%   words, the trace is performed along the two-dimensional slices
%   defined by dimensions I and J. It is possible to implement tensor
%   multiplication as an outer product followed by a contraction.
%
%   Examples
%   X = sptenrand([4 3 2],10); Y = sptenrand([3 2 4],10);
%   Z1 = ttt(X,Y,1,3); %<-- Normal tensor multiplication
%   Z2 = contract(ttt(X,Y),1,6); %<-- Outer product + contract
%   norm(Z1-Z2) %<-- Should be zero
%
%   See also SPTENSOR, SPTENSOR/TTT.
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
if x.size(i) ~= x.size(j)
    error('Must contract along equally sized dimensions');
end

% Error checking
if i == j
    error('Must contract along two different dimensions');
end

% Easy case - returns a scalar
if ndims(x) == 2
    tfidx = (x.subs(:,1) == x.subs(:,2)); % find diagonal entries
    y = sum(x.vals(tfidx));
    return;
end

% Remaining dimensions after contract
remdims = setdiff(1:ndims(x),[i j]);

% Find index of values on diagonal
indx = find(x.subs(:,i) == x.subs(:,j));

% Let the constructor sum up the entries
y = sptensor(x.subs(indx,remdims),x.vals(indx),x.size(remdims));

% Check if result should be dense
if nnz(y) > 0.5 * prod(y.size)
    % Final result is a *dense* tensor
    y = tensor(y);
end
