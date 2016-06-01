% LRR | ROSL | Robust Orthonormal Subspace Learning (Shu et al. 2014)
% process_video('LRR', 'ROSL', 'dataset/demo.avi', 'output/demo_LRR-ROSL.avi');

K = 1; % The initialiation of the subspace dimension
tol = 1e-5;
maxIter = 30;
lambda = 1e-1; %2e-3;
[~,~,E_hat,A_hat] = inexact_alm_rosl(M,K,lambda,tol,maxIter);
L = A_hat;
S = E_hat;