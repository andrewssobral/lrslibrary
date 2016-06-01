% PRMF (Wang et al. 2012)
% process_video('RPCA', 'PRMF', 'dataset/demo.avi', 'output/demo_PRMF.avi');
X = normalize(M);
rk = 2;
lambdaU = 1;
lambdaV = 1;
tol = 1e-2;
[P, Q] = RPMF(X, rk, lambdaU, lambdaV, tol);
L = P * Q;
%S = abs(X - P * Q);
S = X - L;