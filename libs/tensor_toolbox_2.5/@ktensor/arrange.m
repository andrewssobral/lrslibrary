function X = arrange(X,foo)
%ARRANGE Arranges the rank-1 components of a ktensor.
%
%   ARRANGE(X) normalizes the columns of the factor matrices and then sorts
%   the ktensor components by magnitude, greatest to least.
%
%   ARRANGE(X,N) absorbs the weights into the Nth factor matrix instead of
%   lambda. 
%
%   ARRANGE(X,P) rearranges the components of X according to the
%   permutation P. P should be a permutation of 1 to NCOMPOMENTS(X). 
%
%   See also KTENSOR, NCOMPONENTS.
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


%% Just rearrange and return if second argument is a permutation
if exist('foo','var') && (length(foo) > 1)
    X.lambda = X.lambda(foo);
    for i = 1 : ndims(X)
        X.u{i} = X.u{i}(:,foo);
    end   
    return;
end

%% Ensure that matrices are normalized
X = normalize(X);

%% Sort
[X.lambda, idx] = sort(X.lambda, 1, 'descend');
for i = 1 : ndims(X)
    X.u{i} = X.u{i}(:,idx);
end

%% Absorb the weight into one factor, if requested
if exist('foo','var')
    r = length(X.lambda);
    X.u{end} = X.u{end} * spdiags(X.lambda,0,r,r);
    X.lambda = ones(size(X.lambda));
end

