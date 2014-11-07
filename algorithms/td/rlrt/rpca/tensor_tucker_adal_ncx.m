function results = tensor_tucker_adal_ncx( data, params )
% Solve Tucker decomposition
%       min_{X,E} ||E||^2
%       s.t. X_i + E = T
%            rank(X_i) <= k_i

% data.T
% params
% X, V are cell arrays of tensors.
%
% output: 
% results.G = core tensor
% results.U{i} = i-th orthorg factor matrix
%
% Algorithm: ADAL

tic;
T = data.X;     %%% use the 'clean' tensor!!
N = length( size(params.X0) );
X = cell( 1, N );
U = cell( 1, N );
V = cell( 1, N );
for i = 1:N
    X{i} = params.X0;
    V{i} = params.V0{i};
end
E = params.E0;
ks = params.k;

% lambda = params.lambda;
mu = params.mu1;
verbose = params.verbose;
tnorm = norm(T);
rel_err = 1;

for iter = 1:params.max_iter
    % solve X_i's
    for i = 1:N
        [X{i}, junk, U{i}] = tensor_hard_thresh( T+mu*V{i}-E, ks(i), i );
    end
    
    % solve E
    Ep = E;
    D = cell( 1, N );
    for i = 1:N
        D{i} = X{i} - T - mu*V{i};
    end
    E = -ten_sum_all(D) / (N+2*mu);
%     E = shrinkage_t( P, mu/N );
    
    % compute optimality stats
    pres = 0;
    tdiff = cell( 1, N );
    for i = 1:N
        tdiff{i} = X{i} + E - T;
    end
    pres = tensor_array_norm(tdiff);
    
    Ediff = E - Ep;
    pres = pres / (sqrt(N)*tnorm);
    dres = norm(Ediff) / norm(Ep);
    Y = ten_sum_all(X) / N;
    rel_err_p = rel_err;
    rel_err = norm(Y-data.X) / norm(data.X);
    err_chg = abs( rel_err_p - rel_err );
    
    % print
    if verbose
    fprintf('Iter: %d,   pinf: %3.2e,   dinf: %3.2e,    rel_err: %3.2e,     err_chg: %3.2e\n', iter, pres, dres, rel_err, err_chg );
    end
    
    if max(pres, dres) < params.opt_tol && err_chg < params.opt_tol
%     if pres < params.opt_tol
        break;
    end
    
    % update Lagrange multipliers
    for i = 1:N
        V{i} = V{i} - tdiff{i}/mu;
    end
    
%     if rem(iter,10) == 0
%         mu = max(mu*0.9, params.mu_min);
%     end
end

%%% compute the core tensor G
UT = cell( 1, N );
for i = 1:N
    UT{i} = U{i}';
end
results.G = ttm(Y, UT, 1:N);
results.X = Y;
results.E = E;
results.V = V;
results.T = T;
results.U = U;
results.iter = iter;
results.cpu = toc;
results.mu = mu;
% results.lambda = lambda;

end
