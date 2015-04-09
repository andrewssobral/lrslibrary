function [A_hat, E_hat, iter, svp, elapsed] = inexact_alm_rpca_with_blws(D, lambda, tol, maxIter, blk)

% This is an exemplar code illustrating how to use block Lanczos with warm start (BLWS)
% (the BL_SVD function) in a host algorithm.
%
% Oct 2009
% This matlab code implements the inexact augmented Lagrange multiplier 
% method for Robust PCA. Reference: 
%
% Zhouchen Lin, Minming Chen, and Yi Ma, The Augmented Lagrange Multiplier 
% Method for Exact Recovery of Corrupted Low-Rank Matrix, Technical Report 
% UILU-ENG-09-2215, UIUC, October 2009 (arXiv: 1009.5055). 
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
% blk - indicate whether to use BLWS
% 
% Initialize A,E,Y,u
% while ~converged 
%   minimize (inexactly, update A and E only once)
%     L(A,E,Y,u) = |A|_* + lambda * |E|_1 + <Y,D-A-E> + mu/2 * |D-A-E|_F^2;
%   Y = Y + \mu * (D - A - E);
%   \mu = \rho * \mu;
% end
%
% Reference: Zhouchen Lin and Siming Wei, A Block Lanczos with Warm Start 
% Technique for Accelerating Nuclear Norm Minimization Algorithms, 
% arxiv: 1012.0365.
%
% Bug report: zclin2000@hotmail.com
%


%addpath PROPACK;

%elapsed = tic;

[m, n] = size(D);

if nargin < 2
    lambda = 1 / sqrt(max(m,n));
elseif lambda == -1
    lambda = 1 / sqrt(max(m,n));
end

if nargin < 3
    tol = 1e-6;
elseif tol == -1
    tol = 1e-6;
end

if nargin < 4
    maxIter = 1000;
elseif maxIter == -1
    maxIter = 1000;
end

if nargin < 5
    blk = 0;
elseif blk == -1
    blk = 0;
end
    

% initialize
Y = D;
norm_two = lansvd(Y, 1, 'L');
norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;

A_hat = zeros( m, n);
E_hat = zeros( m, n);
mu = 1.25/norm_two; % this one can be tuned
mu_bar = mu * 1e7;
rho = 1.3;         % this one can be tuned
d_norm = norm(D, 'fro');

iter = 0;
total_svd = 0;
converged = false;
stopCriterion = 1;
sv = 10;
d = min(m, n);

mark = false;
Block_mark = false ;

while ~converged       
    iter = iter + 1;
    
    temp_T = D - A_hat + (1/mu)*Y;
    E_hat = max(temp_T - lambda/mu, 0);
    E_hat = E_hat+min(temp_T + lambda/mu, 0);
    
    if Block_mark==true
%         uv_temp=[U_temp; V_temp]/sqrt(2);
        [U S V]=BL_SVD(D - E_hat + (1/mu)*Y, U_temp, V_temp, 1);
    else  
          [U S V] = lansvd(D - E_hat + (1/mu)*Y, sv, 'L');
    end
    
    if mark==true
        U_temp = U;
        V_temp = V;
        Block_mark=true;
    end

    diagS = diag(S);
    svp = length(find(diagS > 1/mu));
    
    len = length(diagS);
    ratio = diagS(1:len-1) ./ diagS(2:len);
    ind = find(ratio > 2);
    if blk && length(ind) > 0
        svp = min(svp, ind(1));
        mark = true;
    end   

    if svp < sv
        sv = min(svp + 1, d);
    else
        sv = min(svp + round(0.05*d), d);
    end
    
    A_hat = U(:, 1:svp) * diag(diagS(1:svp) - 1/mu) * V(:, 1:svp)';    

    total_svd = total_svd + 1;
    
    Z = D - A_hat - E_hat;
    
    Y = Y + mu*Z;
    mu = min(mu*rho, mu_bar);
        
    %% stop Criterion    
    stopCriterion = norm(Z, 'fro') / d_norm;
    if stopCriterion < tol
        converged = true;
    end    
    
    if mod( total_svd, 10) == 0
        disp(['#svd ' num2str(total_svd) ' r(A) ' num2str(svp)...
            ' |E|_0 ' num2str(length(find(abs(E_hat)>0)))...
            ' stopCriterion ' num2str(stopCriterion)]);
    end    
    
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end

%elapsed = toc(elapsed);
