function Y = divide(X,K,epsilon)
%DIVIDE Divide an SPTENSOR by a nonnegative KTENSOR.
%
%  Y = DIVIDE(X,K,EPSILON) divides the sparse tensor X by the 
%  nonnegative ktensor K.  Avoids divide-by-zero errors by dividing 
%  by MIN(EPSILON,K-VALUE) at each nonzero of X.
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

% Assumes K is a nonnegative ktensor

Y = X;

subs = Y.subs;
vals = zeros(size(Y.vals));
R = numel(K.lambda);
N = ndims(Y);
for r = 1:R
    tvals = ones(size(vals)) * K.lambda(r);
    for n = 1:N
        v = K{n}(:,r);
        tvals = tvals .* v(subs(:,n));
    end
    vals = vals + tvals;
end
Y.vals = Y.vals ./ max(epsilon, vals);

return;
