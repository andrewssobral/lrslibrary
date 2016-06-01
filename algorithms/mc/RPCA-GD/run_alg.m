% MC | RPCA-GD | Robust PCA via Gradient Descent (Yi et al. 2016) | 1

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

% subsampling
[numr,numc] = size(M);
I = randi([0 1],numr,numc); % ones(size(M));
Y0 = sparse(M.*I);

[U,V] = rpca_gd(Y0, r, alpha_bnd, params);

L = U*V';
S = M - L;
