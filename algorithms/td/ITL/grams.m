function [Q, R] = grams(A)

% grams  Gram-Schmidt orthogonalization of the columns of A.
% The columns of A are assumed to be linearly independent.
%
% Q = grams(A) returns an m by n matrix Q whose columns are 
% an orthonormal basis for the column space of A.
%
% [Q, R] = grams(A) returns a matrix Q with orthonormal columns
% and an invertible upper triangular matrix R so that A = Q*R.
%
% Warning: For a more stable algorithm, use [Q, R] = qr(A, 0) .

[m, n] = size(A);
Asave = A;
for j = 1:n
  for k = 1:j-1
    mult = (A(:, j)'*A(:, k)) / (A(:, k)'*A(:, k));
    A(:, j) = A(:, j) - mult*A(:, k);
  end
end
for j = 1:n
  if norm(A(:, j)) < sqrt(eps)
    error('Columns of A are linearly dependent.')
  end
  Q(:, j) = A(:, j) / norm(A(:, j));
end
R = Q'*Asave;
