% Online PRMF (Wang et al. 2012)
% process_video('RPCA', 'OPRMF', 'dataset/demo.avi', 'output/demo_OPRMF.avi');
X = normalize(M);
rk = 2;
lambdaU = 1;
lambdaV = 1;
tol = 1e-2;
mask = ones(size(X));
[~, ~, L] = onlineRPMF(X, rk, lambdaU, lambdaV, tol, mask);
S = X - L;