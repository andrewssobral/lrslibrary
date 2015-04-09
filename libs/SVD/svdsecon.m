function [U,S,V] = svdsecon(X,k)
% Input:
% X : m x n matrix
% k : extracts the first k singular values
%
% Output:
% X = U*S*V' approximately (up to k)
%
% Description:
% Does equivalent to svds(X,k) but faster
% Requires that k < min(m,n) where [m,n] = size(X)
% This function is useful if k is much smaller than m and n
% or if X is sparse (see doc eigs)
%
% Vipin Vijayan (2014)

%X = bsxfun(@minus,X,mean(X,2));
[m,n] = size(X);
assert(k <= m && k <= n, 'k needs to be smaller than size(X,1) and size(X,2)');

if  m <= n
    C = X*X';
    [U,D] = eigs(C,k);
    clear C;
    if nargout > 2
        V = X'*U;
        s = sqrt(abs(diag(D)));
        V = bsxfun(@(x,c)x./c, V, s');
        S = diag(s);
    end
else
    C = X'*X; 
    [V,D] = eigs(C,k);
    clear C;
    U = X*V; % convert evecs from X'*X to X*X'. the evals are the same.
    %s = sqrt(sum(U.^2,1))';
    s = sqrt(abs(diag(D)));
    U = bsxfun(@(x,c)x./c, U, s');
    S = diag(s);
end
