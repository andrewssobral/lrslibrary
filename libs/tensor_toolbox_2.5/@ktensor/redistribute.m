function X = redistribute(X,mode)
%REDISTRIBUTE Distribute lambda values to a specified mode. 
%
%   K = REDISTRIBUTE(K,N) absorbs the weights from the lambda vector
%   into mode N. Set the lambda vector to all ones.
%
%   See also KTENSOR, NORMALIZE.
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

for r = 1:length(X.lambda)
    X.u{mode}(:,r) = X.u{mode}(:,r) * X.lambda(r);
    X.lambda(r) = 1;
end
