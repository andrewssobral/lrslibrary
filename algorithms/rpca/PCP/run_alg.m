% PCP (Candes et al. 2009)
% process_video('RPCA', 'PCP', 'dataset/demo.avi', 'output/demo_PCP.avi');
lambda = 1/sqrt(max(size(M))); % default lambda
tol = 1e-5;
[L,S] = PCP(M,lambda,tol);
