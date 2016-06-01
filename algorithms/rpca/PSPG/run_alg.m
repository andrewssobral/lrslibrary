% PSPG (Aybat et al. 2012)
% process_video('RPCA', 'PSPG', 'dataset/demo.avi', 'output/demo_PSPG.avi');
stdev = 1;
tol = 0.05;
L = pspg(M,stdev,tol);
S = M - L;