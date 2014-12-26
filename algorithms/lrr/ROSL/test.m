%close all; 
clear all; clc;

%%
load('dataset/trafficdb/traffic_patches.mat');
V = im2double(imgdb{100});
%show_3dvideo(V);

%%% Matrix-based algorithms
[M, m, n, p] = convert_video3d_to_2d(V);
%show_2dvideo(M,m,n);

%%
K = 1;      %   The initialiation of the subspace dimension
tol = 1e-5; 
maxIter = 100;

%%
lambda = 1e-1; 
[Dict, alpha, E_hat, A_hat] = inexact_alm_rosl(M, K, lambda, tol, maxIter);
%%
X = A_hat; %abs(M - A_hat);
show_2dvideo(X(:,1:10),m,n);


%%
L = 51;     %   The size of submatrix 
lambda = 4e-3; %   this weight paramter can be 0.5~3e-2   (1) Lobby:  8e-3  (2) face: 4e-3
[Dict, alpha, E_hat, A_hat] = inexact_alm_rosl_subsampling(M, K, L, lambda, tol, maxIter);

%%
M_hat = Dict*alpha;
%%
show_2dvideo(A_hat,m,n);

%%
for i = 1:51
  K = i;      %   The initialiation of the subspace dimension
  %L = i;     %   The size of submatrix 
  tol = 1e-5; 
  maxIter = 30;
  lambda = 2e-3; 
  [Dict, alpha, E_hat, A_hat] = inexact_alm_rosl(M, K, lambda, tol, maxIter);
  M_hat = Dict*alpha;
  show_2dvideo(M_hat(:,1:3),m,n);
end