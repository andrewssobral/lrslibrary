function B = qmult(A,method)
%QMULT Pre-multiply matrix by random orthogonal matrix.
%   QMULT(A) returns Q*A where Q is a random real orthogonal matrix
%   from the Haar distribution of dimension the number of rows in A.
%   Special case: if A is a scalar then QMULT(A) is the same as QMULT(EYE(A)).
%   QMULT(A,METHOD) specifies how the computations are carried out.
%   METHOD = 0 is the default, while METHOD = 1 uses a call to QR,
%   which is much faster for large dimensions, even though it uses more flops.

%   Called by RANDCOLU, RANDCORR, RANDJORTH, RANDSVD.

%   Reference:
%   G. W. Stewart, The efficient generation of random
%   orthogonal matrices with an application to condition estimators,
%   SIAM J. Numer. Anal., 17 (1980), 403-409.
%
%   Nicholas J. Higham
%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 1.4.4.2 $  $Date: 2005/11/18 14:15:22 $

if nargin < 2 || isempty(method); method = 0; end

[n, m] = size(A);

%  Handle scalar A.
if max(n,m) == 1
   n = A;
   A = eye(n);
end

if method == 1
   [Q,R] = qr(randn(n));
   B = Q*diag(sign(diag(R)))*A;
   return
end

d = zeros(n,1);

for k = n-1:-1:1

    % Generate random Householder transformation.
    x = randn(n-k+1,1);
    s = norm(x);
    sgn = mysign(x(1));
    s = sgn*s;
    d(k) = -sgn;
    x(1) = x(1) + s;
    beta = s*x(1);

    % Apply the transformation to A.
    y = x'*A(k:n,:);
    A(k:n,:) = A(k:n,:) - x*(y/beta);

end

% Tidy up signs.
for i=1:n-1
    A(i,:) = d(i)*A(i,:);
end
A(n,:) = A(n,:)*mysign(randn);
B = A;
