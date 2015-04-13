function f = flops_lu(n)
% FLOPS_LU   Flops for LU decomposition.
% FLOPS_LU(n) returns the number of flops to compute lu(rand(n,n)).
% The matrix is assumed to be symmetric positive definite, so that no pivoting is required.

% Number of flops for the algorithm in Cormen et al:
% Number of multiplies+adds is:
%   sum(k=1:n) sum(i=k+1:n) sum(j=k+1:n) 3
% = sum(k=1:n) 3(n-k)^2
% = sum(k=0:(n-1)) 3 k^2 
% = n(n-0.5)(n-1)
% Number of divides is: n-1
f = n*(n-0.5)*(n-1) + (n-1)*flops_div;
