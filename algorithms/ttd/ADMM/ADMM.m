%%% ADMM
% A = L + S + Z
% L: low-rank
% S: sparse
% Z: stochastic or deterministic pertubation
function results = ADMM(M)
N = 3;
[mm,nn] = size(M);
g2_max = norm(M(:),inf); % Inf-norm 
g3_max = norm(M); % L2-norm
g2 = 0.15*g2_max;
g3 = 0.15*g3_max;
MAX_ITER = 100;
ABSTOL   = 1e-4;
RELTOL   = 1e-2;
lambda = 1;
rho = 1/lambda;
Z = zeros(mm,nn); S = zeros(mm,nn); L = zeros(mm,nn);
z = zeros(mm,N*nn); U = zeros(mm,nn);
fprintf('\n%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
    'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
for k = 1:MAX_ITER
    B = avg(Z, S, L) - M./N + U;

    % x-update
    Z = (1/(1+lambda))*(Z - B);
    S = prox_l1(S - B, lambda*g2);
    L = prox_matrix(L - B, lambda*g3, @prox_l1);

    % (for termination checks only)
    x = [Z S L];
    zold = z;
    z = x + repmat(-avg(Z, S, L) + M./N, 1, N);

    % u-update
    U = B;
    
    % diagnostics, reporting, termination checks
    results.objval(k)   = objective(Z, g2, S, g3, L);
    results.r_norm(k)   = norm(x - z,'fro');
    results.s_norm(k)   = norm(-rho*(z - zold),'fro');
    results.eps_pri(k)  = sqrt(mm*nn*N)*ABSTOL + RELTOL*max(norm(x,'fro'), norm(-z,'fro'));
    results.eps_dual(k) = sqrt(mm*nn*N)*ABSTOL + RELTOL*sqrt(N)*norm(rho*U,'fro');
    
    if k == 1 || mod(k,1) == 0
      fprintf('%4d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
        results.r_norm(k), results.eps_pri(k), results.s_norm(k), results.eps_dual(k), results.objval(k));
    end
    if results.r_norm(k) < results.eps_pri(k) && results.s_norm(k) < results.eps_dual(k)
       break;
    end
end
results.iter = k;
results.Z = Z;
results.S = S;
results.L = L;
M_hat = L + S + Z;

error = norm(M_hat(:)-M(:))/norm(M(:));
disp(['Error: ' num2str(error)]);
end
