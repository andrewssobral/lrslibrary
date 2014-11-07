function nrm = norm(A)
%NORM Frobenius norm of a ktensor.
%
%   NORM(T) returns the Frobenius norm of a ktensor.
%
%   See also KTENSOR.
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


% Retrieve the factors of A
U = A.u;

% Compute the matrix of correlation coefficients
coefMatrix = A.lambda * A.lambda';
for i = 1:ndims(A)
  coefMatrix = coefMatrix .* (U{i}'*U{i});
end

nrm = sqrt(abs(sum(coefMatrix(:))));

return;
