function [sigma,bnd] = bdsqr(alpha,beta)

% BDSQR: Compute the singular values and bottom element of
%        the left singular vectors of a (k+1) x k lower bidiagonal 
%        matrix with diagonal alpha(1:k) and lower bidiagonal beta(1:k),
%        where length(alpha) = length(beta) = k.
%
% [sigma,bnd] = bdsqr(alpha,beta)
%
% Input parameters:
%   alpha(1:k)   : Diagonal elements.
%   beta(1:k)    : Sub-diagonal elements.
% Output parameters:
%   sigma(1:k)  : Computed eigenvalues.
%   bnd(1:k)    : Bottom elements in left singular vectors.

% Below is a very slow replacement for the BDSQR MEX-file.

% warning('PROPACK:NotUsingMex','Using slow matlab code for bdsqr.')
k = length(alpha);
if min(size(alpha)') ~= 1  | min(size(beta)') ~= 1
  error('alpha and beta must be vectors')
elseif length(beta) ~= k
  error('alpha and beta must have the same lenght')
end    
B = spdiags([alpha(:),beta(:)],[0,-1],k+1,k);
[U,S,V] = svd(full(B),0);
sigma = diag(S);
bnd = U(end,1:k)';


