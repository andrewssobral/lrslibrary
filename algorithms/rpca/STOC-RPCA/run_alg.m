% RPCA | STOC-RPCA | Online Robust PCA via Stochastic Optimization (Feng et al. 2013)
% process_video('RPCA', 'STOC-RPCA', 'dataset/demo.avi', 'output/demo_STOC-RPCA.avi');
lambda1 = 1/sqrt(max(size(M)));
lambda2 = lambda1;
nrank = size(M,2);
[L,S] = stoc_rpca(M, lambda1, lambda2, nrank);
