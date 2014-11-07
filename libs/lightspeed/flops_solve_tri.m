function f = flops_solve_tri(n,m,k)
% FLOPS_SOLVE_TRI   Flops for triangular left division.
% FLOPS_SOLVE_TRI(T,b) returns the number of flops for solve_tri(T,b).
% FLOPS_SOLVE_TRI(n,m,k) returns the number of flops for 
% solve_tril(tril(rand(n,m)),rand(n,k)).
%
% Example:  (n=2,m=2,k=1)
%  [g;h] = [a 0; b c]\[e; f]  has
%  g = e/a
%  h = (f - b*g)/c
%  which is 2 multiply+add and 2 divisions = 18 flops.

if nargin == 2
	T = n;
	b = m;
  f = flops_solve_tri(rows(T),cols(T),cols(b));
  return;
end
if n ~= m
	error('n ~= m case is not implemented');
end
% lower triangular case:
% number of multiplies+adds is
% sum(i=1..n) sum(k=1..i-1) 2 = sum(i=1..n) 2*(i-1) = n^2-n
% number of divides is n
f = (n*n + n*(flops_div-1))*k;
