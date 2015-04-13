function [A_hat, E_hat, iter] = exact_alm_rpca(D, lambda, tol, maxIter)

% Oct 2009
% This matlab code implements the augmented Lagrange multiplier method for
% Robust PCA.
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
%   minimize
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

[m,n] = size(D);

if(nargin < 2) lambda = 1 / sqrt(m); end
if(nargin < 3) tol = 1e-7; elseif(tol == -1) tol = 1e-7; end
if(nargin < 4) maxIter = 1000; elseif(maxIter == -1) maxIter = 1000; end

% initialize
Y = sign(D);
norm_two = norm(Y, 2); %norm_two = lansvd(Y, 1, 'L');
norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;

A_hat = zeros( m, n);
E_hat = zeros( m, n);
dnorm = norm(D, 'fro');
tolProj = 1e-6 * dnorm;
total_svd = 0;
mu = .5/norm_two; % this one can be tuned
rho = 6;          % this one can be tuned

iter = 0;
converged = false;
stopCriterion = 1;
sv = 5;
svp = sv;
while ~converged       
    iter = iter + 1;
    % solve the primal problem by alternative projection
    primal_converged = false;
    primal_iter = 0;
    sv = sv + round(n * 0.1);
    while primal_converged == false
        temp_T = D - A_hat + (1/mu)*Y;
        temp_E = max( temp_T - lambda/mu,0) + min( temp_T + lambda/mu,0); 
        
        %if choosvd(n, sv) == 1
        %    [U S V] = lansvd(D - temp_E + (1/mu)*Y, sv, 'L');
        %else
            %[U,S,V] = svd(D - temp_E + (1/mu)*Y, 'econ');
            [U,S,V] = svdecon(D - temp_E + (1/mu)*Y); % fastest
        %end
        
        diagS = diag(S);
        svp = length(find(diagS > 1/mu));
        if svp < sv
            sv = min(svp + 1, n);
        else
            sv = min(svp + round(0.05*n), n);
        end
        temp_A = U(:,1:svp)*diag(diagS(1:svp)-1/mu)*V(:,1:svp)';    
        
        if norm(A_hat - temp_A, 'fro') < tolProj && norm(E_hat - temp_E, 'fro') < tolProj
            primal_converged = true;
        end
        A_hat = temp_A;
        E_hat = temp_E;
        primal_iter = primal_iter + 1;
        total_svd = total_svd + 1;
    end
        
    Z = D - A_hat - E_hat;        
    Y = Y + mu*Z;
    mu = rho * mu;
    
    %% stop Criterion    
    stopCriterion = norm(Z, 'fro') / dnorm;
    if stopCriterion < tol
        converged = true;
    end    
    
    disp(['Iteration' num2str(iter) ' #svd ' num2str(total_svd) ' r(A) ' num2str(svp)...
        ' |E|_0 ' num2str(length(find(abs(E_hat)>0)))...
        ' stopCriterion ' num2str(stopCriterion)]);
    
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end

if nargin == 5
    fclose(fid);
end

