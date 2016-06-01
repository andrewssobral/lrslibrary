function [U, V] = rpca_gd(Y, r, alpha, params)
% [U, V] = RPCA_GD(Y, r, alpha, params)
% Robust PCA via Non-convex Gradient Descent
%
% Y : A sparse matrix to be decomposed into a low-rank matrix M and a sparse
% matrix S. Unobserved entries are represented as zeros.
% r : Target rank
% alpha : An upper bound of max sparsity over the columns/rows of S
% params : parameters for the algorithm
%   .step_const : Constant for step size (default .5)
%   .max_iter : Maximum number of iterations (default 30)
%   .tol : Desired Frobenius norm error (default 2e-4)
%   .incoh : Incoherence of M (default 5)
%
% Output:
% U, V : M=U*V' is the estimated lowrank matrix
%
% By:
% Xinyang Yi, Dohyung Park, Yudong Chen, Constantine Caramanis
% {yixy,dhpark,constantine}@utexas.edu, yudong.chen@cornell.edu


% Default parameter settings
step_const = .5;
max_iter   = 30;
tol        = 2e-4;
do_project = 0;
gamma      = 1;
incoh      = 5;

% Read paramter settings
if isfield(params,'gamma')      gamma = params.gamma; end
if isfield(params,'incoh')      incoh = params.incoh; end
if isfield(params,'step_const') step_const = params.step_const; end
if isfield(params,'max_iter')   max_iter = params.max_iter; end
if isfield(params,'tol')        tol= params.tol; end
if isfield(params,'do_project') do_project = params.do_project; end

% Library paths
%addpath PROPACK;
%addpath MinMaxSelection;

% Setting up
err  = zeros(1,max_iter);
time = zeros(1,max_iter);
Ynormfro = norm(Y,'fro');
[d1, d2] = size(Y);

is_sparse  = issparse(Y);
if is_sparse
    [I, J, Y_vec] = find(Y);
    n = length(Y_vec);
    obs_ind = sub2ind([d1,d2], I, J);
    col = [0; find(diff(J)); n];
    p = n/d1/d2;
    if p>0.9
        is_sparse = 0;
        Y = full(Y);
    end
else
    p = 1;
end

%% Phase I: Initialization
t1 = tic; t = 1;

% Initial sparse projection
fprintf('Initial sparse projection; time %f \n', toc(t1));
alpha_col = alpha; alpha_row = alpha;
S = Tproj_partial(Y, gamma*p*alpha_col, gamma*p*alpha_row);

% Initial factorization
fprintf('Initial SVD; time %f \n', toc(t1));

[U,Sig,V] = lansvd((Y-S)/p,r,'L');
U = U(:,1:r) * sqrt(Sig(1:r,1:r));
V = V(:,1:r) * sqrt(Sig(1:r,1:r));

% Projection
if do_project
    const1 = sqrt(4*incoh*r/d1)*Sig(1,1);
    const2 = sqrt(4*incoh*r/d2)*Sig(1,1);
    U = U .* repmat(min(ones(d1,1),const1./sqrt(sum(U.^2,2))),1,r);
    V = V .* repmat(min(ones(d2,1),const2./sqrt(sum(V.^2,2))),1,r);
end

% Compute the initial error
err(t)  = inf;

time(t) = toc(t1);

%% Phase II: Gradient Descent
steplength = step_const / Sig(1,1);

if is_sparse
    YminusUV = sparse(I, J, 1, d1, d2, n);
else
    YminusUV = zeros(d1, d2);
end

fprintf('Begin Gradient descent\n');
converged = 0;
while ~converged
    
    fprintf('Iter no. %d err %e time %f \n', t, err(t), time(t));
    t = t + 1;
    
    %%
    if is_sparse
        UVobs_vec = compute_X_Omega(U, V, obs_ind);
        %UVobs_vec = partXY(U', V', I, J, n)';
        YminusUV = sparse(I, J, Y_vec-UVobs_vec, d1, d2, n); clearvars UVobs_vec;
    else
        YminusUV = Y - U*V';
    end
    %err(t) = norm(YminusUV-S, 'fro')/Ynormfro;
    
    %% Sparse Projection for S
    S = Tproj_partial(YminusUV, gamma*p*alpha_col, gamma*p*alpha_row);
    E = YminusUV - S;
    clearvars S;
    
    %% Gradient Descent for U and V
    Unew = U + steplength * (E * V) /p - steplength/16*U*(U'*U-V'*V);
    Vnew = V + steplength * (U' * E)' /p - steplength/16*V*(V'*V-U'*U);
    
    %% Projection
    if do_project
        Unew = Unew .* repmat(min(ones(d1,1),const1./sqrt(sum(Unew.^2,2))),1,r);
        Vnew = Vnew .* repmat(min(ones(d2,1),const2./sqrt(sum(Vnew.^2,2))),1,r);
    end
    
    U = Unew;
    V = Vnew;
    
    %% Compute error
    err(t) = norm(E, 'fro')/Ynormfro;
    time(t) = toc(t1);
    
    %% Convergence check
    if (t >= max_iter)
        converged = 1;
        fprintf('Maximum iterations reached.\n');
    end
    if (err(t) <= max(tol,eps))
        converged = 1;
        fprintf('Target error reached.\n');
    end
    if (err(t) >= err(t-1) - eps)
        converged = 1;
        fprintf('No improvement.\n');
    end
    
end

err  = err(1:t);
time = time(1:t);
