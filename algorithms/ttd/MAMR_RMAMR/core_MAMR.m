function [ A,E,iter1 ] = core_MAMR( D, lambda, I, J, tol, maxiter )
%core_MAMR implements the inexact augmented lagrange multiplier
% method for matrix recovery with erase and sparse noise
%   D    - m*n matrix of observations
%
% lambda - weight on sparse error term in the cost function
%
% I, J   - two dimension index of the regions containing values 
%
% tol    - tolerance for stopping criterion
%    -DEFAULT 1e-3 if omitted or -1
%
% maxiter - maximum number of iterations
%    -DEFAULT 500 if omitted or -1
%
% Model:
%   min |A|_* +lambda*|E|_1
%   subj P(A+E)=P(D);
% Copyright:Xinchen YE, Tianjin University, 2014

[m,n] = size(D);

if nargin < 5
    tol = 1e-3;
elseif tol == -1
    tol = 1e-3;
end
if nargin < 6
    maxiter = 500;
elseif maxiter == -1
    maxiter = 500;
end

rho = 1.2; % iteration setp
norm_two = lansvd(D, 1, 'L'); %computes the 1 largest singular value
muk = 1000/norm_two; % penalty initialization 
d_norm=norm(D,'fro');

Ek=zeros(m,n);Yk=zeros(m,n);
iter1=0;
converged1=false;
sv =1;
while ~converged1
    iter1 = iter1+1;
    
    [U, S, V]=lansvd(D-Ek+(1/muk)*Yk,sv,'L');
%     [U,S,V] = svd(D-Ek+(1/muk)*Yk);
    diagS = diag(S);
    diagS = diagS(1:sv);
    svn = length(find(diagS > 1/muk));
    svp = svn;
    
    ratio = diagS(1:end-1)./diagS(2:end);
    [max_ratio, max_idx] = max(ratio);
    if max_ratio > 2
        svp = min(svn, max_idx);
    end
    if svp < sv %|| iter < 10
        sv = min(svp + 1, n);
    else
        sv = min(svp + 10, n);
    end
    Ak = U(:,1:svp)*diag(diagS(1:svp)-1/muk)*V(:,1:svp)';
%     Ak = U*(shrink(S,1/muk))*V';
    
    Ek = MtOmega(shrink(D-Ak+(1/muk)*Yk,lambda/muk),I,J,m,n)+...
        D-Ak+(1/muk)*Yk-MtOmega(D-Ak+(1/muk)*Yk,I,J,m,n);
    
    Yk = Yk + muk * (D-Ak-Ek);
    muk = rho * muk;
   
     stopCriterion = norm(D-Ak-Ek, 'fro') / d_norm;
      disp([ ' r(F) ' num2str(rank(Ak))...
            ' |E|_0 ' num2str(length(find(abs(Ek)>0)))...
            ' stopCriterion ' num2str(stopCriterion)  ' iter1 ' num2str(iter1) ' mu ' num2str(muk)]);
    if stopCriterion < tol
        converged1 = true;
    end   
    if ~converged1 && iter1 >= maxiter
        disp('Maximum iterations reached');
        converged1 = true;
    end
end

A=Ak;
E=Ek;
end

