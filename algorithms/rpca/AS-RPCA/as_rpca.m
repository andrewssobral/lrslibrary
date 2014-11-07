function [A,E,obj] = as_rpca(D,lambda,k,rho,display)
% This routine solves the following optimization problem,
% min |A|_*+lambda*|E|_1
% s.t., X = Z+E, rank(A) < = k
% reference: Active Subspace: Towards Scalable Low-Rank Learning, Neural
% Computation, 2012.

[d n] = size(D);
if nargin<5
    display = false;
end

if nargin<4
    rho = 1.1;
end

if nargin<3
     k = max(1,round(min(d,n)/10));
end

if nargin<2
    lambda = 1/sqrt(min(d,n));
end

tol = 1e-8*norm(D,'fro');
maxIter = 1000;
max_mu = 1e10;
mu = 1/norm(D,2);
%% Initializing optimization variables
% intialize
J = zeros(k,n);
E = sparse(d,n);

Y = zeros(d,n);
%% Start main loop
iter = 0;
objs=[];
while iter<maxIter
    iter = iter + 1;
    dey = D - E + Y/mu;
    temp = dey*J';
    %update Q
    [U,sigma,V] = svd(temp,'econ');
    Q = U*V';
    
    %update J
    temp = Q'*dey;

    [U,sigma,V] = svd(temp,'econ');
    sigma = diag(sigma);
    svp = length(find(sigma>1/mu));
    if svp>=1
        sigma = sigma(1:svp)-1/mu;
    else
        svp = 1;
        sigma = 0;
    end
    J = U(:,1:svp)*diag(sigma)*V(:,1:svp)';
    
    %update E
    A = Q*J;
    temp = D - A + Y/mu;
    E = max(0,temp - lambda/mu)+min(0,temp + lambda/mu);
    
    leq = D - A - E;
    stopC = norm(leq,'fro');
    if display
        E_norm = sum(sum(abs(D - A)));
        obj = sum(sigma)+lambda*E_norm;
        objs = [objs,obj];
    end
    if display && (iter==1 || mod(iter,10)==0 || stopC<tol)
        disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e') ...
            ',obj=' num2str(obj) ',stopALM=' num2str(stopC,'%2.3e')]);
    end
    if stopC<tol 
        break;
    else
        Y = Y + mu*leq;
        mu = min(max_mu,mu*rho);
    end
end
if display
    figure;
    plot(1:length(objs),objs);
    xlabel('#iteration');
    ylabel('objective function value');
end
