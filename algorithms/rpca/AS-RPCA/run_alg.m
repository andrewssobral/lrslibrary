% AS-RPCA: Active Subspace: Towards Scalable Low-Rank Learning (Liu and Yan, 2012)
% process_video('RPCA', 'AS-RPCA', 'dataset/demo.avi', 'output/demo_AS-RPCA.avi');
lambda = 1/sqrt(min(size(M)));
[L,S] = as_rpca(M,lambda);
