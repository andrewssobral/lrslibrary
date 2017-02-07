%{
load('dataset/trafficdb/traffic_patches.mat');
[M,m,n,p] = convert_video3d_to_2d(im2double(imgdb{100}));
out = run_algorithm('MC', 'RPCA-GD', M, [])
show_results(M.*out.Omega,out.L,out.S,out.O,p,m,n);
%}

% Parameter setting
% d1 = 5e4; d2 = 6e4; % matrix size
% r  = 5;  % rank
% p  = 0.005;  % observation probability
% alpha = 0.05; % sparsity of S* (expected)
% Another parameter setting with full observation
r = 1;
alpha = 0.1;
% Decomposition via Gradient Descent
% algorithm paramters
params.step_const = 0.5; % step size parameter for gradient descent
params.max_iter   = 30;  % max number of iterations
params.tol        = 2e-4;% stop when ||Y-UV'-S||_F/||Y||_F < tol
% alpha_bnd is some safe upper bound on alpha, 
% that is, the fraction of nonzeros in each row of S (can be tuned)
gamma = 1.5;
alpha_bnd = gamma*alpha;
Y0 = sparse(M.*Omega); % subsampling
[U,V] = rpca_gd(Y0, r, alpha_bnd, params);
L = U*V'; % low-rank
S = M - L; % sparse
