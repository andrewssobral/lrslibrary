%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = pdct_operator(picks,perm)

% Define A*x and A'*y for a partial DCT matrix A
% Input:
%            n = interger > 0
%        picks = sub-vector of a permutation of 1:n
% Output:
%        A = struct of 2 fields
%            1) A.times: A*x
%            2) A.trans: A'*y

if exist('dct','file')
    DCT = @dct;  IDCT = @idct;
elseif exist('dct2','file')
    DCT = @dct2; IDCT = @idct2;
else
    error('DCT functions not found');
end

A.times = @(x) pdct_n2m( DCT,x,picks,perm);
A.trans = @(y) pdct_m2n(IDCT,y,picks,perm);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = pdct_n2m(DCT,x,picks,perm)

% Calculate y = A*x,
% where A is m x n, and consists of m rows of the 
% n by n discrete-cosine transform (DCT) matrix
% with columns permuted by perm.
% The row indices are stored in picks.

x = x(:);
tx = DCT(x(perm));
y = tx(picks);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = pdct_m2n(IDCT,y,picks,perm)

% Calculate x = A'*y,
% where A is m x n, and consists of m rows of the 
% n by n inverse discrete-cosine transform (IDCT)
% matrix with columns permuted by perm.
% The row indices are stored in picks.

n = length(perm);
tx = zeros(n,1);
tx(picks) = y;
x = zeros(n,1);
x(perm) = IDCT(tx);
