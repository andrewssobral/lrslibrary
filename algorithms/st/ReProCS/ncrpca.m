function [L_t, S_t, iters, frob_err] = ncrpca(M, true_r, EPS, MAX_ITER, EPS_S, incoh, TOL)

% This matlab code implements Non-convex Robust PCA (NcRPCA)
% Input:
% M = given low rank+sparse matrix to be decomposed
% true_r = maximum rank of the low rank rank component
% EPS (optional) = convergence threshold for ||M-(L_t+S_t)||_F; default is 1e-3
% MAX_ITER (optional) = maximum iterations for NcRPCA; default is 51
% EPS_S (optional) = threshold for removing small entries in the sparse component; default is 1e-3
% incoh (optional) = incoherence of the low rank component; default is 1
% TOL (optional) = tolerance for relative error in ||M-(L_t+S_t)||_F in consecutive iterations; default is 1e-1
% Output:
% M_t = thresholded M at each iteration
% L_t = rank-k approximation of M_t
% S_t = sparse component, computed as M-M_t
% iters = number of iteration of NcRPCA
% frob_err = ||M-(L_t+S_t)||_F at each iteration

if nargin < 7, TOL = 1e-1; end
if nargin < 6, incoh = 1; end
if nargin < 5, EPS_S = 1e-3; end
if nargin < 4, MAX_ITER = 51; end
if nargin < 3, EPS = 1e-3; end

%addpath code_ncrpca/PROPACK;
%addpath PROPACK;
frob_err(1) = inf;
[~, n] = size(M);
t = 1;
idx = [];
thresh_const = incoh; % threshold constant: can be tuned depending on incoherence
thresh_red = 0.9; % parameter to reduce the threshold constant: can be tuned
r_hat = 1; % initial rank for stagewise algorithm
L_t = zeros(size(M));
Sig_t = lansvd(M,1,'L');
D_t = M-L_t;
thresh = thresh_const*Sig_t/sqrt(n);
idx = unique([find(abs(D_t) > thresh); idx]);
S_t = zeros(size(M));
S_t(idx) = D_t(idx); % initial thresholding
if max(idx(:))==0
    idx = [];
end
while frob_err(t)/norm(M, 'fro')>=EPS && t<MAX_ITER % convergence check
    if ~mod(t, 10)
        %fprintf('Iter no. %d\n', t);
    end
    t = t+1;
    [U_t, Sig_t, V_t] = lansvd(M-S_t, r_hat+1, 'L');
    %[U_t, Sig_t, V_t] = svds(M-S_t, r_hat+1); % use this if not using propack
    L_t= U_t(:,1:r_hat)*Sig_t(1:r_hat,1:r_hat)*V_t(:,1:r_hat)';
    D_t = M-L_t;
    thresh = (thresh_const/sqrt(n))*Sig_t(r_hat+1, r_hat+1); % use n instead of sqrt(n) for thresholding less aggresively
    idx = unique([find(abs(D_t) > thresh); idx]);
    S_t(idx) = D_t(idx);
    frob_err(t) = norm(M-(L_t+S_t), 'fro');
    if ((frob_err(t-1)-frob_err(t))/frob_err(t-1) <= TOL) && r_hat<true_r
%         r_hat = r_hat+1; % use this for incrementally updating rank by 1
        sig_t = lansvd(M-S_t, true_r, 'L'); % svd function from propack
        ratio_sig = sig_t(r_hat+1:end)./[sig_t(r_hat+2:end); sig_t(end)];
        [~, mx_idx] = max(ratio_sig);
        r_hat = r_hat+mx_idx; % update rank for the next stage
    elseif ((frob_err(t-1)-frob_err(t))/frob_err(t-1) <= TOL) && r_hat==true_r
        thresh_const = thresh_const*thresh_red; % tune threshold
    end
end
S_t(abs(S_t)<EPS_S) = 0; % threshold to remove small errors and obtain sparse component
iters = length(frob_err)-1; % no. of iters. of ncrpca = length(frob_err)-1
end
