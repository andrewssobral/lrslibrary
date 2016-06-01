% RPCA | RegL1-ALM | Low-Rank Matrix Approximation under Robust L1-Norm (Zheng et al. 2012)
% process_video('RPCA', 'RegL1-ALM', 'dataset/demo.avi', 'output/demo_RegL1-ALM.avi');
W = ones(size(M));
r = 1;
lambda = 1e-3;
rho = 1.2;
maxIterIN = 1;
signM = 0;
%[M_est,U_est,V_est,L1_error] = ...
L = RobustApproximation_M_UV_TraceNormReg(M,W,r,lambda,rho,maxIterIN,signM);
S = M - L;
