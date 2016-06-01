% =========================================================================
% 
% Testing Robust PCA on synthetic data
% Using a Non-convex Gradient Descent Algorithm
% 
% =========================================================================

% Parameter setting
% d1 = 5e4; d2 = 6e4; % matrix size
% r  = 5;  % rank
% p  = 0.005;  % observation probability
% alpha = 0.05; % sparsity of S* (expected)

% Another parameter setting with full observation
d1 = 5e3; d2 = 6e3; 
r  = 5; 
p  = 1;
alpha = 0.1;


%% Generate a problem instance: Y0 = P_Omega( U0*V0' + S0 )
fprintf('Generating Y0 = U0*V0'' + S0\n');
ncol = floor(p*d1);
n  = d2*ncol;
U0 = randn(d1,r)/sqrt(d1);
V0 = randn(d2,r)/sqrt(d2);
Y0 = sparse([],[],[],d1,d2,n);
I0 = zeros(n,1);
J0 = zeros(n,1);
X0 = zeros(n,1);
for j=1:d2

    I0(((j-1)*ncol+1):j*ncol) = sort(randsample(d1,ncol));
    J0(((j-1)*ncol+1):j*ncol) = j;
    X0(((j-1)*ncol+1):j*ncol) = U0(I0(((j-1)*ncol+1):j*ncol),:)*V0(j,:)';
    
    IS   = randsample((j-1)*ncol+1:j*ncol,floor(alpha*ncol));
    X0(IS) = (r/sqrt(d1*d2)) * (rand(numel(IS),1)*5) .* sign(randn(numel(IS),1));

end
Y0 = sparse(I0,J0,X0,d1,d2);
fprintf('%d x %d, %d observed, %d corrupted, rank %d\n', d1, d2, n, floor(alpha*n), r);


%% Decomposition via Gradient Descent
% algorithm paramters
params.step_const = 0.5; % step size parameter for gradient descent
params.max_iter   = 30;  % max number of iterations
params.tol        = 2e-4;% stop when ||Y-UV'-S||_F/||Y||_F < tol

% alpha_bnd is some safe upper bound on alpha, 
% that is, the fraction of nonzeros in each row of S (can be tuned)
gamma = 1.5;
alpha_bnd = gamma*alpha;   

t1 = tic;
[U, V] = rpca_gd(Y0, r, alpha_bnd, params);
t2 = toc(t1);

%% Subspace errors
fprintf('\n');
fprintf('Computation time: %f\n', t2);
fprintf('Subspace error: %f\n', max(subspace(U, U0), subspace(V, V0)) );