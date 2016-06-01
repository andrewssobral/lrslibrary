% APG Partial (Lin et al. 2009)
% process_video('RPCA', 'APG_PARTIAL', 'dataset/demo.avi', 'output/demo_APG_PARTIAL.avi');
lambda = 1/sqrt(max(size(M))); % default lambda
[L,S] = partial_proximal_gradient_rpca(M,lambda);
