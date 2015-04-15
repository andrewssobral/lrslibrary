clear;clc;

% Data generation
m = 100;
n = 100;
r = 0.05*m;

U = randn(m, r);
V = randn(r, n);
M = U*V;

% Mixture noise
ind = randperm(m*n);
num1 = 0.1*m*n;
num2 = 0.2*m*n;
num3 = 0.7*m*n;

N1 = zeros(m, n);
N1(ind(1:num1)) = 50*rand(1, num1)-25;
Y = M + N1;
N2 = 1*randn(1, num2);
N3 = 0.1*randn(1, num3);
Y(ind(num1+1:num1+num2)) = Y(ind(num1+1:num1+num2)) + N2;
Y(ind(num1+num2+1:end)) = Y(ind(num1+num2+1:end)) + N3;

% Perform MoG-RPCA
param.mog_k = 3;
param.lr_init = 'SVD';
param.maxiter = 100;
param.initial_rank = 2*r;
param.tol = 1e-3;
    
lr_prior.a0 = 1e-6;
lr_prior.b0 = 1e-6;

mog_prior.mu0 = 0;
mog_prior.c0 = 1e-3;
mog_prior.d0 = 1e-3;
mog_prior.alpha0 = 1e-3;
mog_prior.beta0 = 1e-3;

[lr_model, mog_model, r] = mog_rpca(Y, param, lr_prior, mog_prior);

L = lr_model.U*lr_model.V';
rre = norm(L-M,'fro')/norm(M,'fro');

display(['Relative reconstruction error: ', num2str(rre)]);
display(['Estimated rank: ', num2str(r)]);