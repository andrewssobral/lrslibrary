function [lr_model, mog_model, r] = mog_rpca(Y, param, lr_prior, mog_prior)

% MoG-RPCA
%   
% Inputs:
%    Y          ----  input data matrix
%    param      ----  input parameters
%       param.maxiter      : number of iterations allowed (default: 100)
%       param.tol          : stop criterion (default: 1e-4)
%       param.mog_k        : number of Gaussians in the noise component (default: 3)
%       param.lr_init      : method for initializing the low-rank component
%                             'SVD'  : using SVD (default)
%                             'rand' : random initialization
%       param.initial_rank : initial rank of the low-rank component (default: full rank)
%    lr_prior   ----  hyperparameters of the low-rank component
%    mog_prior  ----  hyperparameters of the MoG noise component
% Outputs:
%    lr_model   ---- estimated model parameters of the low-rank component
%    mog_model  ---- estimated model parameters of the MoG noise component
%    r          ---- estimated rank of the low-rank matrix
%
% Written by Qian Zhao (if you have any questions/comments/suggestions, please contact me: timmy.zhaoqian@gmail.com)
%
% Code is free for research purposes, but please contact us if for other purposes. In addition, I require that you do not laugh at my code :)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Qian Zhao, Deyu Meng, Zongben Xu, Wangmeng Zuo, Lei Zhang. Robust Principal Component Analysis with Complex Noise. ICML, 2014."
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters initialization
[m,n] = size(Y);
mn = m*n;

if (~isfield(param,'maxiter'))
    maxiter = 100;
else
    maxiter = param.maxiter;
end

if (~isfield(param,'tol'))
    tol = 1e-4;
else
    tol = param.tol;
end

if (~isfield(param,'mog_k'))
    mog_k = 3;
else
    mog_k = param.mog_k;
end

if (~isfield(param,'lr_init'))
    lr_init = 'SVD';
else
    lr_init = param.lr_init;
end

if (~isfield(param,'initial_rank'))
    initial_rank = min([m,n]);
else
    initial_rank = param.initial_rank;
end

clear param;

k = mog_k;

% Low-rank model hyperparameters
if nargin < 3
    lr_prior.a0 = 1e-6;
    lr_prior.b0 = 1e-6;
end

% MoG model hyperparameters
if nargin < 4
    mog_prior.mu0 = 0;
    mog_prior.c0 = 1e-6;
    mog_prior.d0 = 1e-6;
    mog_prior.alpha0 = 1e-6;
    mog_prior.beta0 = 1e-6;
end

% Low-rank model initialization
Y2sum = sum(Y(:).^2);
scale2 = Y2sum / (mn);
scale = sqrt(scale2);
if strcmp(lr_init, 'SVD')   % SVD initialization
    %[u, s, v] = svd(Y, 'econ');
    [u, s, v] = svdecon(Y); % fastest
    r = initial_rank;      
    U = u(:,1:r)*(s(1:r,1:r)).^(0.5);
    V = (s(1:r,1:r)).^(0.5)*v(:,1:r)'; 
    V = V';     
    Sigma_U = repmat( scale*eye(r,r), [1 1 m] );
    Sigma_V = repmat( scale*eye(r,r), [1 1 n] );
    gammas = scale*ones(r,1);
elseif strcmp(lr_init, 'rand')  % Random initialization
    r = initial_rank;  
    U = randn(m,r) * sqrt(scale);
    V = randn(n,r) * sqrt(scale);
    Sigma_U = repmat( scale*eye(r,r), [1 1 m] );
    Sigma_V = repmat( scale*eye(r,r), [1 1 n] );
    gammas = scale*ones(r,1);
end
lr_model.U = U;
lr_model.V = V;
lr_model.Sigma_U = Sigma_U;
lr_model.Sigma_V = Sigma_V;
lr_model.gammas = gammas;
L = U*V';
    
% MoG model Initialization
E = Y - L;
mog_model.R = R_initialization(E(:)', k);
nk = sum(mog_model.R,1);
nxbar = E(:)'*mog_model.R;
mog_model.c = mog_prior.c0 + nk/2;
mog_model.beta = mog_prior.beta0 + nk;
mog_model.d = mog_prior.d0 + 0.5*((E(:)'.^2*mog_model.R)+mog_prior.beta0*mog_prior.mu0^2-1./mog_model.beta.*(nxbar+mog_prior.beta0*mog_prior.mu0).^2);
mog_model.R = reshape(mog_model.R, m, n, k);
mog_model.mu = 1./mog_model.beta.*(mog_prior.beta0.*mog_prior.mu0+nxbar);

% Main loop
for iter=1:maxiter
    L_old = L;
    
    % LR update
    [lr_model, r, E_YminusUV, E_YminusUV_2] = lr_update(Y, lr_model, mog_model, r, lr_prior);
    L = lr_model.U*lr_model.V';
    
    % mog update
    mog_model = mog_vmax(mog_model, mog_prior, E_YminusUV, E_YminusUV_2);
    mog_model = mog_vexp(mog_model, E_YminusUV, E_YminusUV_2);
    
    % Convergence check
    if norm(L-L_old,'fro')/norm(L_old,'fro') < tol
        break;
    end  
end
[~, label] = max(mog_model.R, [], 3);
mog_model.label = label;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lr_model, r, E_YminusUV, E_YminusUV_2] = lr_update(Y, lr_model, mog_model, r, lr_prior)
[m, n] = size(Y);

a0 = lr_prior.a0;
b0 = lr_prior.b0;

U = lr_model.U;
V = lr_model.V;
Sigma_U = lr_model.Sigma_U;
Sigma_V = lr_model.Sigma_V ;
gammas = lr_model.gammas;

R = mog_model.R;
c = mog_model.c;
d = mog_model.d;
mu = mog_model.mu;

k = length(mu);

tau = c./d;
Gam = diag(gammas);
Rtau = reshape(reshape(R,m*n,k)*tau',m,n);
Rtaumu = reshape(reshape(R,m*n,k)*(tau.*mu)',m,n);
RtauYmu = Rtau.*Y - Rtaumu;

% Update U
re_Sigma_V = reshape(Sigma_V, r*r, n);
diagsU = zeros(r,1);
temp_U = zeros(r,r,m);
for i=1:m
    Sigma_U(:,:,i) = ( reshape( re_Sigma_V*Rtau(i,:)', r, r ) + bsxfun(@times, V', Rtau(i,:))*V + Gam )^(-1);
    U(i,:) = (RtauYmu(i,:)*V) * Sigma_U(:,:,i);
    diagsU = diagsU + diag( Sigma_U(:,:,i) );
    temp_U(:,:,i) = Sigma_U(:,:,i)+U(i,:)'*U(i,:);
end

% Update V
re_Sigma_U = reshape(Sigma_U, r*r, m);
diagsV = zeros(r,1);
temp_V = zeros(r,r,n);
for j=1:n
    Sigma_V(:,:,j) = ( reshape( re_Sigma_U*Rtau(:,j), r, r ) + bsxfun(@times, U', Rtau(:,j)')*U + Gam )^(-1);
    V(j,:) = (RtauYmu(:,j)'*U) * Sigma_V(:,:,j);
    diagsV = diagsV + diag( Sigma_V(:,:,j) );
    temp_V(:,:,j) = Sigma_V(:,:,j)+V(j,:)'*V(j,:);
end

% Update gammas
gammas = ( 2*a0 + m + n )./( 2*b0 + sum(U.*U)'+ diagsU + sum(V.*V)' + diagsV);

% Prune redundant dimesions
dim_thr = 1e2;
max_gamma = min(gammas) * dim_thr;
if sum(find(gammas > max_gamma)),
    indices = find(gammas <= max_gamma);    
    U = U(:,indices);
    V = V(:,indices);
    gammas = gammas(indices);
    Sigma_U = Sigma_U(indices,indices,:);
    Sigma_V = Sigma_V(indices,indices,:);
    temp_U = temp_U(indices,indices,:);
    temp_V = temp_V(indices,indices,:);
    r = size(U,2);
end

lr_model.U = U;
lr_model.V = V;
lr_model.Sigma_U = Sigma_U;
lr_model.Sigma_V = Sigma_V;
lr_model.gammas = gammas;

E_YminusUV = Y - U*V';
E_YminusUV_2 = Y.^2 - 2.*Y.*(U*V') + reshape(reshape(temp_U,r*r,m)'*reshape(temp_V,r*r,n),m,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mog_model = mog_vmax(mog_model, mog_prior, E_YminusUV, E_YminusUV_2)
alpha0 = mog_prior.alpha0;
beta0 = mog_prior.beta0;
mu0 = mog_prior.mu0;
c0 = mog_prior.c0;
d0 = mog_prior.d0;
R = mog_model.R;

[m, n] = size(E_YminusUV);
k = size(R,3);

nxbar = reshape(E_YminusUV, 1, m*n)*reshape(R, m*n, k);

nk = sum(reshape(R, m*n, k),1);
alpha = alpha0+nk;
beta = beta0+nk;
c = c0+nk/2;
mu = bsxfun(@times,bsxfun(@plus,beta0*mu0,nxbar),1./beta);
d = d0 + 0.5*( reshape(E_YminusUV_2, 1, m*n)*reshape(R, m*n, k) + beta0*mu0^2 -1./beta.*(nxbar+beta0*mu0).^2 );

mog_model.alpha = alpha;
mog_model.beta = beta;
mog_model.mu = mu;
mog_model.c = c;
mog_model.d = d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  mog_model = mog_vexp(mog_model, E_YminusUV, E_YminusUV_2)
alpha = mog_model.alpha;
beta = mog_model.beta;
mu = mog_model.mu;
c = mog_model.c;
d = mog_model.d;

[m, n] = size(E_YminusUV);
k = length(mu);
Ex = reshape(E_YminusUV, m*n, 1);
Ex2 = reshape(E_YminusUV_2, m*n, 1);

tau = c./d;
EQ = zeros(m*n, k);
for i=1:k
    EQ(:,i) = 1/beta(i) + tau(i)*mu(i)^2 + tau(i)*Ex2 - 2*tau(i)*mu(i)*Ex;
end

Elogtau = psi(0, c) - log(d);
Elogpi = psi(0, alpha) - psi(0, sum(alpha));

logRho = (bsxfun(@minus,EQ,2*Elogpi+Elogtau-log(2*pi)))/(-2);
logR = bsxfun(@minus,logRho,logsumexp(logRho,2));
R = exp(logR);

mog_model.logR = reshape(logR, m, n, k);
mog_model.R = reshape(R, m, n, k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = R_initialization(X, k)
n = size(X, 2);
idx = randsample(n,k);
m = X(:,idx);
[~,label] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1);
[u,~,label] = unique(label);
while k ~= length(u)
    idx = randsample(n,k);
    m = X(:,idx);
    [~,label] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1);
    [u,~,label] = unique(label);
end
R = full(sparse(1:n,label,1,n,k,n));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%