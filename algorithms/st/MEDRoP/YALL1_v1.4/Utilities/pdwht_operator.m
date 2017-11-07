%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = pdwht_operator(picks,perm)

% Define A*x and A'*y for a partial DWHT matrix A
% Input:
%            n = interger of a power of 2
%        picks = sub-vector of a permutation of 1:n
% Output:
%        A = struct of 2 fields
%            1) A.times: A*x
%            2) A.trans: A'*y

A.times = @(x) pfwht_n2m(x,picks,perm);
A.trans = @(y) pfwht_m2n(y,picks,perm);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = pfwht_n2m(x,picks,perm)

% Calculate y = A*x,
% where A is m x n, and consists of m rows of the 
% n by n discrete-Walsh-Hadamard transform matrix
% with permuted columns by perm. 
% The row indices are stored in picks.

x = x(:);
n = length(x);
tx = fWHtrans(x(perm))*sqrt(n);
y = tx(picks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = pfwht_m2n(y,picks,perm)

% Calculate x = A'*y,
% where A is m x n, and consists of m rows of the 
% n by n discrete-Walsh-Hadamard transform matrix
% with permuted columns by perm.
% The row indices are stored in picks.

n = length(perm);
tx = zeros(n,1);
tx(picks) = y/sqrt(n);
x = zeros(n,1);
x(perm) = ifWHtrans(tx);
