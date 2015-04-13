function f = flops_inv(n)
% FLOPS_INV    Flops for matrix inversion.
% FLOPS_INV(n) returns the number of flops for inv(rand(n,n)).
% For n=1,2 the number of flops is exact.
% For n>2, the number of flops is estimated by assuming the matrix is 
% decomposed via LU or Cholesky and then inverted by back-substitution.
% See flops_solve.

% Written by Tom Minka

if n == 2
  % inv([a b; c d]) = [d -b; -c a]/(ad-bc)
  % 3 flops for determinant, 1 divide, 4 multiplies, 1 negation
  % thanks to Zhang Xiaoying for pointing out this special case.
  f = 3+flops_div + 5;
else
  f = flops_solve(n,n,n);
end


