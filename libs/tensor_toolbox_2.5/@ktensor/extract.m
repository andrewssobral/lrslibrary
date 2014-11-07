function new_X = extract(X,idx)
%EXTRACT Creates a new ktensor with only the specified components.
%
%   Y = EXTRACT(X,S) selected the subset of components in X as defined by
%   S. It should be the case that S is a subset of [1,...,NCOMPONENTS(X)].
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


%% Set-up
N = ndims(X);
%% Extract
new_lambda = X.lambda(idx);
new_U = cell(N,1);
for i = 1 : N
    new_U{i} = X.u{i}(:,idx);
end
new_X = ktensor(new_lambda, new_U);

