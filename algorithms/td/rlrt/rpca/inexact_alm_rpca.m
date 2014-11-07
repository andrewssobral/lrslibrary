function [A_hat E_hat iter] = inexact_alm_rpca(D, lambda, tol, maxIter, mu)

% Oct 2009
% This matlab code implements the inexact augmented Lagrange multiplier 
% method for Robust PCA.
%
% D - m x n matrix of observations/data (required input)
%
% lambda - weight on sparse error term in the cost function
%
% tol - tolerance for stopping criterion.
%     - DEFAULT 1e-7 if omitted or -1.
%
% maxIter - maximum number of iterations
%         - DEFAULT 1000, if omitted or -1.
% 
% Initialize A,E,Y,u
% while ~converged 
%   minimize (inexactly, update A and E only once)
%     L(A,E,Y,u) = |A|_* + lambda * |E|_1 + <Y,D-A-E> + mu/2 * |D-A-E|_F^2;
%   Y = Y + \mu * (D - A - E);
%   \mu = \rho * \mu;
% end
%
% Minming Chen, October 2009. Questions? v-minmch@microsoft.com ; 
% Arvind Ganesh (abalasu2@illinois.edu)
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing

% addpath PROPACK;
global use_propack;

[m n] = size(D);

if nargin < 2
    lambda = 1 / sqrt(m);
end

if nargin < 3
    tol = 1e-3;%1e-7;
elseif tol == -1
    tol = 1e-7;
end

if nargin < 4
    maxIter = 1000;
elseif maxIter == -1
    maxIter = 1000;
end

% initialize
Y = D;
norm_two = lansvd(Y, 1, 'L');
norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;

A_hat = zeros( m, n);
E_hat = zeros( m, n);
% mu = 1.25/norm_two; % this one can be tuned
mu_bar = mu * 1e7;
rho = 1.5;          % this one can be tuned
d_norm = norm(D, 'fro');

iter = 0;
total_svd = 0;
converged = false;
stopCriterion = 1;
d = min(m,n);
sv = min(10, d);
while ~converged       
    iter = iter + 1;
    
    temp_T = D - A_hat + (1/mu)*Y;
    Ep = E_hat;
    E_hat = max(temp_T - lambda/mu, 0);
    E_hat = E_hat+min(temp_T + lambda/mu, 0);

    if choosvd(n, sv) == 1 && use_propack
        [U S V] = lansvd(D - E_hat + (1/mu)*Y, sv, 'L');
    else
        [U S V] = svd(D - E_hat + (1/mu)*Y, 'econ');
    end
    diagS = diag(S);
    svp = length(find(diagS > 1/mu));
    if svp < sv
        sv = min(svp + 1, d);
    else
        sv = min(svp + round(0.05*d), d);
    end
    
    A_hat = U(:, 1:svp) * diag(diagS(1:svp) - 1/mu) * V(:, 1:svp)';    

    total_svd = total_svd + 1;
    
    Z = D - A_hat - E_hat;
    Ediff = E_hat - Ep;
    
    Y = Y + mu*Z;
%     mu = min(mu*rho, mu_bar);
        
    %% stop Criterion    
    stopCriterion = norm(Z, 'fro') / d_norm;
    dres = norm(Ediff,'fro') / norm(Ep,'fro');
    if max(stopCriterion,dres) < tol && iter > 10
        converged = true;
    end    
    
%     if mod( total_svd, 10) == 0
%         disp(['#svd ' num2str(total_svd) ' r(A) ' num2str(rank(A_hat))...
%             ' |E|_0 ' num2str(length(find(abs(E_hat)>0)))...
%             ' stopCriterion ' num2str(stopCriterion)]);
%     end    
    
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end
