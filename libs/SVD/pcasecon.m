function [U,T,mu] = pcasecon(X,k)
% Input:
% X : m x n matrix
% Each column of X is a feature vector
%
% Output:
% X = U*T approximately (up to k)
%
% Description:
% Principal Component Analysis (PCA) while trying to conserve memory and be faster
% Requires that k <= min(m,n) where [m,n] = size(X)
% This function is useful if k is much smaller than m and n
% or if X is sparse (see doc eigs)
%
% Vipin Vijayan (2014)

mu = mean(X,2);
X = bsxfun(@minus,X,mu);

[m,n] = size(X);
assert(k <= m && k <= n, 'k needs to be smaller than size(X,1) and size(X,2)');

if  m <= n
    C = X*X';
    [U,~] = eigs(C,k);
    clear C;

    T = U'*X;
else
    C = X'*X; 
    [V,D] = eigs(C,k);
    clear C;
    
    U = X*V; % convert evecs from X'*X to X*X'. the evals are the same.
    %s = sqrt(sum(U.^2,1))';
    s = sqrt(abs(diag(D)));
    U = bsxfun(@(x,c)x./c, U, s');
    %S = diag(s);
    %T = S*V';
    T = bsxfun(@(c,vt)c.*vt,s,V');
end

