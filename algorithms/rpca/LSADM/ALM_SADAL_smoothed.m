function out = ALM_SADAL_smoothed(D, opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test ALM in the paper ''Fast Alternating Linearization Methods for
%       Minimizing the Sum of Two Convex Functions'', Donald Goldfarb,
%       Shiqian Ma and Katya Scheinberg, Tech. Report, Columbia University,
%       2009 - 2010. 
%
% Author: Shiqian Ma
% Date  : Apr. 20, 2010 
% IEOR, Columbia University, Copyright (2010)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n] = size(D);
mu = opts.mu; sigma = opts.sigma; rho = opts.rho;
Y = zeros(m,n); gradgY = zeros(m,n);
X = D - Y; Dnorm = norm(D,'fro');
Lambda = zeros(m,n); sv = opts.sv;

for itr = 1: opts.maxitr
    if choosvd(n, sv) == 1
        [U, gamma, V] = lansvd(mu*Lambda-Y+D,sv,'L');
    else
        [U, gamma, V] = svd(mu*Lambda-Y+D, 'econ');
    end
    gamma = diag(gamma);
    gamma_new = gamma-mu*gamma./max(gamma,mu+sigma);
    svp = length(find(gamma > mu));
    if svp < sv
        sv = min(svp + 1, n);
    else
        sv = min(svp + round(0.05*n), n);
    end

    Xp = X;
    X = U*diag(gamma_new)*V';
    Lambda = Lambda - (X+Y-D)/mu;
    muY = mu;
    B = Lambda - (X-D)/muY;
    Yp = Y;

    Y = muY*B - muY*min(rho, max(-rho,muY*B/(sigma+muY)));
    Lambda = Lambda - (X+Y-D)/muY;

    StopCrit = norm(D-X-Y,'fro')/Dnorm;
    relX = norm(X-opts.Xs,'fro')/norm(opts.Xs,'fro');
    relY = norm(Y-opts.Ys,'fro')/norm(opts.Ys,'fro');
%     fprintf('iter: %d, mu: %3.2e, rank(X):%d, relX: %3.2e, relY: %3.2e,crit:%3.2e\n', ...
%         itr, mu, length(find(gamma_new>1e-3)), relX, relY, norm(D-X-Y,'fro')/norm(D,'fro'));
    %fprintf('iter: %d, mu: %3.2e, crit:%3.2e\n', itr, mu, norm(D-X-Y,'fro')/norm(D,'fro'));

    mu = max(opts.muf, mu*opts.eta_mu);
    sigma = max(opts.sigmaf, sigma*opts.eta_sigma);

    if StopCrit < opts.epsilon
        out.X = X; out.Y = Y; out.iter = itr; out.relX = relX; out.relY = relY; out.StopCrit = StopCrit;
        return;
    end
end

out.X = X;
out.Y = Y; 
out.iter = itr; 
out.relX = relX; 
out.relY = relY; 
out.StopCrit = StopCrit;
