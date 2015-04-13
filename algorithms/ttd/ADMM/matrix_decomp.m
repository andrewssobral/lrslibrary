function h = matrix_decomp

%% Problem instance settings

m = 20;
n = 50;
use_cvx = 1; % set to 0 for larger instances

%% Problem data

s = RandStream.create('mt19937ar','seed',5489);
RandStream.setDefaultStream(s);

N = 3;
r = 4;

L = randn(m,r) * randn(r,n);    % low rank
S = sprandn(m,n,0.05);          % sparse
S(S ~= 0) = 20*binornd(1,0.5,nnz(S),1)-10;
V = 0.01*randn(m,n);            % noise

A = S + L + V;

g2_max = norm(A(:),inf);
g3_max = norm(A);
g2 = 0.15*g2_max;
g3 = 0.15*g3_max;

%% CVX

if use_cvx
    tic;

    cvx_begin
        cvx_precision low
        variables X_1(m,n) X_2(m,n) X_3(m,n)
        minimize(0.5*square_pos(norm(X_1,'fro')) + g2*norm(X_2(:),1) + g3*norm_nuc(X_3))
        subject to
            X_1 + X_2 + X_3 == A;
    cvx_end

    h.cvx_toc = toc;
    h.p_cvx = cvx_optval;
    h.X1_cvx = X_1;
    h.X2_cvx = X_2;
    h.X3_cvx = X_3;

    X_2(abs(X_2) < 1e-4) = 0;
    rhat = sum(svd(X_3) > 1e-4);
    fprintf('CVX (vs true):\n');
    fprintf('|V| = %.2f;  |X_1| = %.2f\n', norm(V, 'fro'), norm(X_1,'fro'));
    fprintf('nnz(S) = %d; nnz(X_2) = %d\n', nnz(S), nnz(X_2));
    fprintf('rank(L) = %d; rank(X_3) = %d\n', rank(L), rhat);
end

%% ADMM

MAX_ITER = 100;
ABSTOL   = 1e-4;
RELTOL   = 1e-2;

tic;

lambda = 1;
rho = 1/lambda;

X_1 = zeros(m,n);
X_2 = zeros(m,n);
X_3 = zeros(m,n);
z   = zeros(m,N*n);
U   = zeros(m,n);

fprintf('\n%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
    'r norm', 'eps pri', 's norm', 'eps dual', 'objective');

for k = 1:MAX_ITER

    B = avg(X_1, X_2, X_3) - A./N + U;

    % x-update
    X_1 = (1/(1+lambda))*(X_1 - B);
    X_2 = prox_l1(X_2 - B, lambda*g2);
    X_3 = prox_matrix(X_3 - B, lambda*g3, @prox_l1);

    % (for termination checks only)
    x = [X_1 X_2 X_3];
    zold = z;
    z = x + repmat(-avg(X_1, X_2, X_3) + A./N, 1, N);

    % u-update
    U = B;

    % diagnostics, reporting, termination checks
    h.objval(k)   = objective(X_1, g2, X_2, g3, X_3);
    h.r_norm(k)   = norm(x - z,'fro');
    h.s_norm(k)   = norm(-rho*(z - zold),'fro');
    h.eps_pri(k)  = sqrt(m*n*N)*ABSTOL + RELTOL*max(norm(x,'fro'), norm(-z,'fro'));
    h.eps_dual(k) = sqrt(m*n*N)*ABSTOL + RELTOL*sqrt(N)*norm(rho*U,'fro');
    
    if k == 1 || mod(k,10) == 0
        fprintf('%4d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
            h.r_norm(k), h.eps_pri(k), h.s_norm(k), h.eps_dual(k), h.objval(k));
    end

    if h.r_norm(k) < h.eps_pri(k) && h.s_norm(k) < h.eps_dual(k)
         break;
    end

end

h.admm_toc = toc;
h.admm_iter = k;
h.X1_admm = X_1;
h.X2_admm = X_2;
h.X3_admm = X_3;

fprintf('\nADMM (vs true):\n');
fprintf('|V| = %.2f;  |X_1| = %.2f\n', norm(V, 'fro'), norm(X_1,'fro'));
fprintf('nnz(S) = %d; nnz(X_2) = %d\n', nnz(S), nnz(X_2));
fprintf('rank(L) = %d; rank(X_3) = %d\n', rank(L), rank(X_3));

if use_cvx
    fprintf('\nADMM vs CVX solutions (in Frobenius norm):\n');
    fprintf('X_1: %.2e; X_2: %.2e; X_3: %.2e\n', ...
        norm(h.X1_cvx - X_1,'fro'), norm(h.X2_cvx - X_2,'fro'), norm(h.X3_cvx - X_3,'fro'));
end

end

function x = avg(varargin)
    N = length(varargin);
    x = 0;
    for k = 1:N
        x = x + varargin{k};
    end
    x = x/N;
end

function p = objective(X_1, g_2, X_2, g_3, X_3)
    p = norm(X_1,'fro').^2 + g_2*norm(X_2(:),1) + g_3*norm(svd(X_3),1);
end
