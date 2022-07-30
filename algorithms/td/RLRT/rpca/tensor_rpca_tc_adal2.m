function results = tensor_rpca_tc_adal2( data, params )
% Solve HO-RPCA with missing data
%       min_{X,E} \sum_i ||X_(i)||_* + \lambda1*||E||_1
%       s.t. A_c(X + E) = T_c
% 
% Reformulate as
%       min \sum_i ||Xi,(i)||_* + \lambda1*||E||_1
%       s.t.    A_c(Xi+E) = T_c, for i = 1,...,N.
%
% data.b = observations (vector)
% params
% X, V are cell arrays of tensors.
%
% Algorithm: I-ADAL

tic;
b = data.b;
bnorm = norm(b);
Omega = data.linInd;
N = length( size(params.X0) );
X = cell( 1, N );
U = cell( 1, N );
V = cell( 1, N );
for i = 1:N
    X{i} = params.X0;
    V{i} = zeros( size(b) );
end
E = params.E0;

lambda = params.lambda;
% lambdaS = [ 1/5, 1/5, 1/5, 1 ];
lambdaS = ones( 1, N );
mu = params.mu1;
eta = 1/1.1;
rel_err = 1;

% assign operators
Aprod = @(X_param)A_select( X_param, Omega );
Atprod = @(b_param)At_select( b_param, Omega, size(E) );

for iter = 1:params.max_iter
    % solve X_i's
    for i = 1:N
        grad = Atprod( Aprod(X{i}+E) - b - mu*V{i} );
        P = X{i} - eta*grad;
        [X{i}, junk, U{i}] = tensor_shrinkage( P, mu*eta*lambdaS(i), i );
    end
    
    % solve E
    Ep = E;
    D = cell( 1, N );
    for i = 1:N
        D{i} = Aprod(X{i}) - b - mu*V{i};
    end
    P = Atprod( -ten_sum_all(D)/N );
    E = shrinkage_t( P, mu*lambda/N );
    
    % compute optimality stats
    pres = 0;
    tdiff = cell( 1, N );
    for i = 1:N
        tdiff{i} = Aprod(X{i} + E) - b;
        pres = pres + norm(tdiff{i})^2;
    end
    pres = sqrt(pres) / (sqrt(N)*bnorm);
    
    Ediff = E - Ep;
    dres = norm(Aprod(Ediff)) / norm(Aprod(Ep));
    
    Y = ten_sum_all(X) / N;
    rel_err_p = rel_err;
    rel_err = norm(Y-data.X) / norm(data.X);
    imp = abs(rel_err-rel_err_p) / rel_err_p;
    
    % print
    if params.verbose
        fprintf('Iter: %d,   pinf: %3.2e,   dinf: %3.2e,    rel_err: %3.2e, imp: %3.2e\n', iter, pres, dres, rel_err, imp );
    end
    
    if max(pres, dres) < params.opt_tol %|| imp < 1e-4
%     if pres < params.opt_tol
        break;
    end
    
    % update Lagrange multipliers
    for i = 1:length(V)
        V{i} = V{i} - tdiff{i}/mu;
    end
    
end

results.X = Y;
results.E = E;
results.V = V;
results.b = b;
results.U = U;
results.iter = iter;
results.cpu = toc;
results.mu = mu;
results.lambda1 = lambda;
results.lambdaS = lambdaS;

end

function AX = A_select( X, Omega )
% X could be a tensor
% AX is a vector

X = double(X);
AX = X(Omega);

end

function X = At_select( b, Omega, sz )

X = zeros( sz );
X(Omega) = b;
X = tensor(X);
end