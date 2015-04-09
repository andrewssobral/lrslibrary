function [mu, Lambda] = sampleHyper(U, m, invW_0, beta_0, nu_0, r)
    U_bar = mean(U)'; 
    Sigma_bar = cov(U); 
    % Sample Lambda_u, W is simplified due to hyperparameter settings
    W = inv(invW_0 + m * Sigma_bar + beta_0 * m / (beta_0 + m) * (U_bar * U_bar'));
    W = (W + W')/ 2;
    Lambda =  (wishrnd(W, nu_0 + m));
    % Sample mu_u
    tempmu = (m * U_bar) / (beta_0 + m);
    lam = chol(inv((beta_0 + m) * Lambda));
    lam = lam';
    mu = lam * randn(r, 1) + tempmu;
end