% Dual RPCA (Lin et al. 2009)
% process_video('RPCA', 'DUAL', 'dataset/demo.avi', 'output/demo_DUAL.avi');
lambda = 1/sqrt(max(size(M))); % default lambda
[L,S] = dual_rpca_2(M,lambda);
