% RPCA | MoG-RPCA | Mixture of Gaussians RPCA (Zhao et al. 2014)
% process_video('RPCA', 'MoG-RPCA', 'dataset/demo.avi', 'output/demo_MoG-RPCA.avi');

%{
clear, clc;
load('dataset/trafficdb/traffic_patches.mat');
V = im2double(imgdb{100});
[M,m,n,p] = convert_video3d_to_2d(V);
%}

r = 1;
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

[lr_model, mog_model, r] = mog_rpca(M, param, lr_prior, mog_prior);

L = lr_model.U*lr_model.V';
S = M - L;

%{
show_2dvideo(M,m,n);
show_2dvideo(L,m,n);
show_2dvideo(S,m,n);
show_2dvideo(hard_threshold(S),m,n);
%}