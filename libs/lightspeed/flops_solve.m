function f = flops_solve(n,m,c)
% FLOPS_SOLVE    Flops for matrix left division.
% FLOPS_SOLVE(a,b) returns the number of flops for a\b.
% FLOPS_SOLVE(n,m,c) returns the number of flops for rand(n,m)\rand(m,c).

if nargin == 2
	a = n;
	b = m;
  f = flops_solve(rows(a),cols(a),cols(b));
  return;
end
if n == m
  if n == 1
    % scalar division
    f = c*flops_div;
  elseif 0
		% invert using cholesky (see inv_posdef)
		f = flops_chol(n) + 2*flops_solve_tri(n,m,c);
	else
		% invert using LU decomposition
		% L has unit diagonal so n divisions are avoided when back-substituting
		f = flops_lu(n) + 2*flops_solve_tri(n,m,c) - n*flops_div;
  end
elseif n > m
  % this comes from Ax=b, x = (A'*A)\(A'*b)
  f = flops_mul(m,n,m) + flops_mul(m,n,c) + flops_solve(m,m,c);
else
  % this comes from Ax=b, x = A'*(A*A')\b
  f = flops_mul(n,m,n) + flops_mul(m,n,c) + flops_solve(n,n,c);
end
