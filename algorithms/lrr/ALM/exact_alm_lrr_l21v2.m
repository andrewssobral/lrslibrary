function [Z_hat, E_hat] = exact_alm_lrr_l21v2(D, A, lambda, tol, maxIter,display)

% Aug 2013
% This matlab code implements the Exact ALM algorithm for
% min_{Z,E}  |Z|_* + lambda |E|_2,1  s.t.  D = AZ + E
%
% D - m x n matrix of observations/data (required input)
% A - m x k matrix of the dictionary (required input) 

% lambda - weight on sparse error term in the cost function
%
% tol - tolerance for stopping criterion.
%     - DEFAULT 1e-7 if omitted or -1.
%
% maxIter - maximum number of iterations
%         - DEFAULT 1000, if omitted or -1.
% 
[m n] = size(D);
k = size(A,2);


if nargin < 4 || isempty(tol)
    tol = 1e-7;
end

if nargin < 5 || isempty(maxIter)
    maxIter = 1000;
end

if nargin<6 || isempty(display)
    display = false;
end

maxIter_primal = 10000;
% initialize
Y = sign(D);
norm_two = norm(Y,2);
norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;

W = zeros(k,n);

Z_hat = zeros(k,n);
E_hat = zeros(m,n);
%parameters
dnorm = norm(D, 'fro');
tolProj1 = 1e-6 * dnorm;

anorm = norm(A,2);
tolProj2 = 1e-6 * dnorm/anorm;

mu = .5/norm_two; % this one can be tuned
rho = 6;          % this one can be tuned

%pre-computation
if m>=k
    inv_ata = inv(eye(k) + A'*A);
else
    inv_ata = eye(k) - A'/(eye(m)+A*A')*A;
end

iter = 0;
while iter < maxIter       
    iter = iter + 1;
    
    % solve the primal problem by alternative projection
    primal_iter = 0;
    
    while primal_iter < maxIter_primal
        primal_iter = primal_iter + 1;
        temp_Z = Z_hat;
        temp_E = E_hat;
        
        %update J
        temp = temp_Z + W/mu;
        [U S V] = svd(temp, 'econ');
       
        diagS = diag(S);
        svp = length(find(diagS > 1/mu));
        diagS = max(0,diagS - 1/mu);
        
        if svp < 0.5 %svp = 0
            svp = 1;
        end
        J_hat = U(:,1:svp)*diag(diagS(1:svp))*V(:,1:svp)';  
        
        % update Z
        temp = J_hat + A'*(D - temp_E) + (A'*Y-W)/mu;
        Z_hat = inv_ata*temp;
        
        %update E
        temp = D - A*Z_hat + Y/mu;
        E_hat =  solve_l1l2(temp, lambda/mu);
        
        if norm(E_hat - temp_E, 'fro') < tolProj1 && norm(Z_hat - temp_Z)<tolProj2
            break;
        end
    end
        
    H1 = D - A*Z_hat - E_hat;
    H2 = Z_hat - J_hat;
    Y = Y + mu*H1;
    W = W + mu*H2;
    mu = rho * mu;
    
    %% stop Criterion    
    stopCriterion = max(norm(H1, 'fro')/dnorm, norm(H2,'fro')/dnorm*anorm);
    if display
        disp(['LRR: Iteration' num2str(iter) '(' num2str(primal_iter) '), mu ' num2str(mu) ', |E|_2,0 ' num2str(sum(sum(E_hat.^2,1)>0))...
        ', stopCriterion ' num2str(stopCriterion)]);
    end
    
    if stopCriterion < tol
        break;
    end    
    
   
end

end



