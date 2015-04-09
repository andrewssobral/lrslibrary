function [Z,E] = adm_lrr(X,lambda,rho,DEBUG)

%------------------------------
% min |Z|_*+lambda*|E|_2,1
% s.t., X = XZ+E
%--------------------------------
% inputs:
%        X -- D*N data matrix, D is the data dimension, and N is the number
%             of data vectors.
% outputs:
%        relChg --- relative changes
%        recErr --- reconstruction errors
if (~exist('DEBUG','var'))
    DEBUG = 0;
end
if nargin < 3
    rho = 1.1;
end
if nargin<2
    lambda = 1;
end
normfX = norm(X,'fro');
tol1 = 1e-4;
tol2 = 1e-5;
maxIter = 1000;
[d n] = size(X);
max_mu = 1e30;
mu = 1e-6;
xtx = X'*X;
inv_x = inv(xtx+eye(n));

opt.tol = tol2;%precision for computing the partial SVD
opt.p0 = ones(n,1);
%% Initializing optimization variables
% intialize
J = zeros(n,n);
Z = zeros(n,n);
E = sparse(d,n);

Y1 = zeros(d,n);
Y2 = zeros(n,n);

sv = 5;
svp = sv;

%% Start main loop
convergenced = 0;
iter = 0;

if DEBUG
    disp(['initial,rank=' num2str(rank(Z))]);
end

while iter<maxIter
    iter = iter + 1;
  
    Em = E;
    Zm = Z;
    
    xmaz = X-X*Z;
    temp = X-X*Z+Y1/mu;
    E = solve_l1l2(temp,lambda/mu);
    
    temp = Z + Y2/mu;
    %[U,sigma,V] = svd(temp,'econ');
    [U,sigma,V] = lansvd(temp,n,n,sv,'L',opt);
    sigma = diag(sigma);
    svp = length(find(sigma>1/mu));
    
    if svp < sv
        sv = min(svp + 1, n);
    else
        sv = min(svp + round(0.05*n), n);
    end
    
    if svp>=1
        sigma = sigma(1:svp)-1/mu;
    else
        svp = 1;
        sigma = 0;
    end
    J = U(:,1:svp)*diag(sigma)*V(:,1:svp)';
    
    Z = inv_x*(xtx-X'*E+J+(X'*Y1-Y2)/mu);
    
    leq1 = xmaz-E;
    leq2 = Z-J;
    relChgZ = norm(Z-Zm,'fro')/normfX;
    relChgE = norm(E-Em,'fro')/normfX;
    relChg = max(relChgE,relChgZ);
    
    recErr = norm(leq1,'fro')/normfX;
    
    convergenced = recErr < tol1 && relChg < tol2;
    
    if DEBUG
        if iter==1 || mod(iter,50)==0 || convergenced
            disp(['iter ' num2str(iter) ',mu=' num2str(mu) ...
                ',rank=' num2str(svp)...
                ',relChg=' num2str(relChg) ' recErr =' num2str(recErr)]);
        end
    end
    if convergenced
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