% NSA v2 (Aybat et al. 2011)
% process_video('RPCA', 'NSA2', 'dataset/demo.avi', 'output/demo_NSA2.avi');
stdev = 1;
tol = 5e-6; % optimality tolerance for stopping_type 1
[L,S] = nsa_v2(M,stdev,tol,1);