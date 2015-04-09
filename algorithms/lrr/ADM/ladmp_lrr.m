function [Z,E] = ladmp_lrr(X,lambda,rho,DEBUG)
% This matlab code implements linearized ADM method for LRR problem
%------------------------------
% min |Z|_*+lambda*|E|_2,1
% s.t., X = XZ+E
%--------------------------------
% inputs:
%        X -- D*N data matrix
% outputs:
%        Z -- N*N representation matrix
%        E -- D*N sparse error matrix
%        relChgs --- relative changes
%        recErrs --- reconstruction errors
%
% created by Risheng Liu on 05/02/2011, rsliu0705@gmail.com
%

%clear global;
clearvars -global M;
global M;% M=Z_k+X'*(X-X*Z_k-E_{k+1}+Y/mu_k)/eta.

%addpath PROPACK;

if (~exist('DEBUG','var'))
    DEBUG = 0;
end
if nargin < 3
    rho = 1.9;
end
if nargin < 2
    lambda = 0.1;
end

normfX = norm(X,'fro');
tol1 = 1e-4;%threshold for the error in constraint
tol2 = 1e-5;%threshold for the change in the solutions
[d n] = size(X);
opt.tol = tol2;%precision for computing the partial SVD
opt.p0 = ones(n,1);

maxIter = 1000;

max_mu = 1e10;
norm2X = norm(X,2);
% mu = 1e2*tol2;
mu = min(d,n)*tol2;

eta = norm2X*norm2X*1.02;%eta needs to be larger than ||X||_2^2, but need not be too large.

%% Initializing optimization variables
% intialize
E = sparse(d,n);
Y = zeros(d,n);
Z = zeros(n, n);

XZ = zeros(d, n);%XZ = X*Z;

sv = 5;
svp = sv;


%% Start main loop
convergenced = 0;
iter = 0;

if DEBUG
    disp(['initial,rank(Z)=' num2str(rank(Z))]);
end

while iter<maxIter
    iter = iter + 1;
    
    %copy E and Z to compute the change in the solutions
    Ek = E;
    Zk = Z;
    
    E = solve_l1l2(X - XZ + Y/mu,lambda/mu);
    
    %-----------Using PROPACK--------------%
    M = Z + X'*(X - XZ - E + Y/mu)/eta;
    
    [U, S, V] = lansvd(M, n, n, sv, 'L', opt);
    %[U, S, V] = lansvd(M, n, n, sv, 'L');
    %[U, S, V] = svd(M,'econ');
      
    S = diag(S);
    svp = length(find(S>1/(mu*eta)));
    if svp < sv
        sv = min(svp + 1, n);
    else
        sv = min(svp + round(0.05*n), n);
    end
    
    if svp>=1
        S = S(1:svp)-1/(mu*eta);
    else
        svp = 1;
        S = 0;
    end

    A.U = U(:, 1:svp);
    A.s = S;
    A.V = V(:, 1:svp);
    
    Z = A.U*diag(A.s)*A.V';

    diffZ = norm(Zk - Z,'fro');
    
    relChgZ = diffZ/normfX;
    relChgE = norm(E - Ek,'fro')/normfX;
    relChg = max(relChgZ,relChgE);

    XZ = X*Z;%introducing XZ to avoid computing X*Z multiple times, which has O(n^3) complexity.
    dY = X - XZ - E;
    recErr = norm(dY,'fro')/normfX;
    
    convergenced = recErr <tol1 && relChg < tol2;
    
    if DEBUG
        if iter==1 || mod(iter,50)==0 || convergenced
            disp(['iter ' num2str(iter) ',mu=' num2str(mu) ...
                ',rank(Z)=' num2str(svp) ',relChg=' num2str(max(relChgZ,relChgE))...
                ',recErr=' num2str(recErr)]);
        end
    end
    if convergenced
%    if recErr <tol1 & mu*max(relChgZ,relChgE) < tol2 %this is the correct
%    stopping criteria. 
        break;
    else
        Y = Y + mu*dY;
        
        if mu*relChg < tol2
            mu = min(max_mu, mu*rho);
        end
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