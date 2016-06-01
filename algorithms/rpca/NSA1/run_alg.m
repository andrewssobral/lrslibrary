% NSA v1 (Aybat et al. 2011)
% process_video('RPCA', 'NSA1', 'dataset/demo.avi', 'output/demo_NSA1.avi');
stdev = 1;
tol = 5e-6; % optimality tolerance for stopping_type 1
L = nsa_v1(M,stdev,tol,1);
S = M - L;