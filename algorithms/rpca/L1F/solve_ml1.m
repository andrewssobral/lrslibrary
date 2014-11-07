% solving 
% min_{E,alpha} |E|_1
% s.t A*alpha+E=D
% where E, D are n*s. A is n*m matrix, alpha is m*s .

function [E, alpha] = solve_ml1(D, A, pinvA)

[n,m] = size(A);
s= size(D, 2);

if nargin<3
   % This step can be improve in future;
   pinvA = pinv(A);
end

% Initialization
alpha = zeros(m,s);
E = zeros(n,s);
Y = zeros(n,s);

tol = 1e-7;
maxIter = 1e6;
rho = 1.1;
max_mu = 1e30;
mu = 1e-6; 

%% Start main loop
iter = 0;
while iter<maxIter
    iter = iter + 1;     
    
    % update E;
    temp_T = D - A*alpha + (1/mu)*Y;
    E = max(temp_T - 1/mu, 0);
    E = E+min(temp_T + 1/mu, 0);
    
    % update alpha;
    temp_T = D + (1/mu)*Y - E;
    alpha = pinvA * temp_T;
    
    leq = D - A*alpha -E;
    stopC = max(max(abs(leq)));
    
    if iter==1 || mod(iter,50)==0 || stopC<tol
        disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e')...
            ',stopALM=' num2str(stopC,'%2.3e')]);
    end
    if stopC<tol 
        break;
    else
    % update Y;
    Y = Y + mu * leq;
    
    % update mu;
    mu =  min(max_mu,mu*rho);
    end
end


