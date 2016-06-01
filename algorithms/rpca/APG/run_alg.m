% APG (Lin et al. 2009)
% process_video('RPCA', 'APG', 'dataset/demo.avi', 'output/demo_APG.avi');
lambda = 1/sqrt(max(size(M))); % default lambda
[L,S] = proximal_gradient_rpca(M,lambda);
