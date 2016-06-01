% SVT (Cai et al. 2008)
% process_video('RPCA', 'SVT', 'dataset/demo.avi', 'output/demo_SVT.avi');
lambda = 1/sqrt(max(size(M))); % default lambda
[L,S] = singular_value_rpca(M,lambda,1e4,0.9,'svd');
