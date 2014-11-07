function [Z,E] = lrr(X,lambda)
% This routine solves the following nuclear-norm optimization problem 
% by using inexact Augmented Lagrange Multiplier, which has been also presented 
% in the paper entitled "Robust Subspace Segmentation 
% by Low-Rank Representation".
%------------------------------
% min |Z|_*+lambda*|E|_2,1
% s.t., X = XZ+E
%--------------------------------
% inputs:
%        X -- D*N data matrix, D is the data dimension, and N is the number
%             of data vectors.
if nargin<2
    lambda = 1;
end
tol = 1e-8;
maxIter = 1e6;
[d n] = size(X);
rho = 1.1;
max_mu = 1e30;
mu = 1e-6;
xtx = X'*X;
inv_x = inv(xtx+eye(n));
%% Initializing optimization variables
% intialize
J = zeros(n,n);
Z = zeros(n,n);
E = sparse(d,n);

Y1 = zeros(d,n);
Y2 = zeros(n,n);
%% Start main loop
iter = 0;
disp(['initial,rank=' num2str(rank(Z))]);
while iter<maxIter
    iter = iter + 1;
    
    temp = Z + Y2/mu;
    [U,sigma,V] = svd(temp,'econ');
 %   [U,sigma,V] = lansvd(temp,30,'L');
    sigma = diag(sigma);
    svp = length(find(sigma>1/mu));
    if svp>=1
        sigma = sigma(1:svp)-1/mu;
    else
        svp = 1;
        sigma = 0;
    end
    J = U(:,1:svp)*diag(sigma)*V(:,1:svp)';
    
    Z = inv_x*(xtx-X'*E+J+(X'*Y1-Y2)/mu);
    
    xmaz = X-X*Z;
    temp = X-X*Z+Y1/mu;
    E = solve_l1l2(temp,lambda/mu);
    
    leq1 = xmaz-E;
    leq2 = Z-J;
    stopC = max(max(max(abs(leq1))),max(max(abs(leq2))));
    if iter==1 || mod(iter,50)==0 || stopC<tol
        disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e') ...
            ',rank=' num2str(rank(Z,1e-3*norm(Z,2))) ',stopALM=' num2str(stopC,'%2.3e')]);
    end
    if stopC<tol 
        break;
    else
        Y1 = Y1 + mu*leq1;
        Y2 = Y2 + mu*leq2;
        mu = min(max_mu,mu*rho);
    end
end

function [E] = solve_l1l2(W,lambda)
n = size(W,2);
E = W;
for i=1:n
    E(:,i) = solve_l2(W(:,i),lambda);
end


function [x] = solve_l2(w,lambda)
% min lambda |x|_2 + |x-w|_2^2
nw = norm(w);
if nw>lambda
    x = (nw-lambda)*w/nw;
else
    x = zeros(length(w),1);
end