function y = logdet(A)
%LOGDET  Logarithm of determinant for positive-definite matrix
% logdet(A) returns log(det(A)) where A is positive-definite.
% This is faster and more stable than using log(det(A)).
% Note that logdet does not check if A is positive-definite.
% If A is not positive-definite, the result will not be the same as log(det(A)).

% Written by Tom Minka
% (c) Microsoft Corporation. All rights reserved.

U = chol(A);
y = 2*sum(log(diag(U)));
