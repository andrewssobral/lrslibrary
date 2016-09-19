function [U_t, SV_t] = ncrmc(M, I)
  %, true_r, p, MAX_ITER, run_time, EPS, EPS_S, incoh, TOL)

%addpath PROPACK
%addpath Mex

TOL = 1e-1;
incoh = 1;
EPS_S = 1e-3;
EPS = 1e-3;
run_time = 100;
MAX_ITER = 70;
true_r = 5;
p = 0.5;

tic;

[m1, m2] = size(M);
[r, c] = find(I);
D_t = M(r + (m1 * (c - 1)));

n = sqrt(m1 * m2);

frob_err(1) = inf;
t = 1;
thresh_const = incoh * true_r; % threshold constant: can be tuned depending on incoherence
thresh_red = 0.9; % parameter to reduce the threshold constant: can be tuned
r_hat = 1; % initial rank for stagewise algorithm

B = D_t;
U_t = zeros(m1, 1);
Sig_t = zeros(1, 1);
V_t = zeros(m2, 1);
[U_t, Sig_t, V_t] = fastsvd(U_t, Sig_t, V_t, (1 / p) * D_t, r, c, 1);
SV_t = Sig_t * V_t';

thresh = thresh_const * Sig_t / n;

S_t = [];

while frob_err(t) >= EPS && t < MAX_ITER  % convergence check

    t = t + 1;

    % thresholding noise
    spL_t = sMatProd(U_t, SV_t, [r, c], 1);
    D_t = B - spL_t;
    idx_s = abs(D_t) >= thresh;
    S_t = D_t;
    S_t(~idx_s) = 0;

    % svp step
    [U_t, Sig_t, V_t] = fastsvd(U_t, Sig_t, V_t, (1 / p) * (D_t - S_t), r, c, r_hat + 1);
    sigma_t = Sig_t(r_hat + 1, r_hat + 1);

    U_t = U_t(:, 1:r_hat);
    Sig_t = Sig_t(1:r_hat, 1:r_hat);
    V_t = V_t(:, 1:r_hat);
    SV_t = Sig_t * V_t';

    % updating threshold
    thresh = (thresh_const / n) * sigma_t; % use n instead of sqrt(n) for thresholding less aggresively

    spL_t = sMatProd(U_t, SV_t, [r, c], 1);
    frob_err(t) = norm(B - (spL_t + S_t), 'fro');

    if ((frob_err(t - 1) - frob_err(t)) / frob_err(t - 1) <= TOL) && r_hat < true_r
        r_hat = r_hat + 1; % use this for incrementally updating rank by 1
    elseif ((frob_err(t - 1) - frob_err(t)) / frob_err(t - 1) <= TOL) && r_hat == true_r
        thresh_const = thresh_const * thresh_red; % tune threshold
    end

end

end
