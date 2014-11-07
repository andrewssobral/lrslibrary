function results = tensor_rpca_tc_adal_ncx( data, params )
% Solve HO-RPCA with missing data
%       min_{X,E} ||E||_1
%       s.t. A_c(X + E) = T_c
%            rank(X_i) <= k_i
% 
% Reformulate as
%       min ||E||_1
%       s.t.    Xi = Y
%               A_c(Y+E) = T_c
%               rank(Xi) <= k_i
%
% data.b = observations (vector)
% params
% X, V are cell arrays of tensors.
%
% Algorithm: ADAL

tic;
b = data.b;
bnorm = norm(b);
Omega = data.linInd;
N = length( size(params.X0) );
X = cell( 1, N );
U = cell( 1, N );
V = cell( 1, N+1 );
for i = 1:N
    X{i} = params.X0;
    V{i} = params.V0{i};
end
V{N+1} = zeros( size(b) );
Y = params.X0;
E = params.E0;
ks = params.k;

mu = params.mu1;
rel_err = 1;

% assign operators
Aprod = @(X_param)A_select( X_param, Omega );
Atprod = @(b_param)At_select( b_param, Omega, size(Y) );

for iter = 1:params.max_iter
    % solve X_i's
    for i = 1:N
        [X{i}, junk, U{i}] = tensor_hard_thresh( Y+mu*V{i}, ks(i), i );
    end
    
    % solve E
    P = Atprod( b - Y(Omega) + mu*V{N+1} );
    E = shrinkage_t( P, mu );
    
    % solve Y
    Yprev = Y;
    tensum = ten_sum_all( X )/mu - ten_sum_all( V(1:N) );
    Y = tensum * mu / N;
    Y(Omega) = (tensum(Omega) + V{N+1} + (b-E(Omega))/mu) / (N/mu + 1/mu);
    
    % compute optimality stats
    pres = 0;
    tdiff = cell( 1, N+2 );
    for i = 1:N
        tdiff{i} = X{i} - Y;
        pres = pres + norm( tdiff{i} )^2;
    end
    tdiff{N+1} = Aprod(Y+E)-b;
    pres = pres + norm(tdiff{N+1})^2;     pres = sqrt(pres);
    Ynorm2 = norm(Y)^2;
    Ydiff = Y - Yprev;
    dres = N*Ydiff;  dres(Omega) = (N-1)*Ydiff(Omega);  dres = norm(dres);
    denom = N*Yprev;    denom(Omega) = (N-1)*Yprev(Omega);     denom = norm(denom);
    pres = pres / sqrt( N*Ynorm2 + bnorm );%norm(b-Aprod(Y))^2 );    % / ynorm;
    dres = dres / denom;
    
    rel_err_p = rel_err;
    rel_err = norm(Y-data.X) / norm(data.X);
    imp = abs(rel_err-rel_err_p) / rel_err_p;
    
    % print
    if params.verbose
        fprintf('Iter: %d,   pinf: %3.2e,   dinf: %3.2e,    rel_err: %3.2e, imp: %3.2e\n', iter, pres, dres, rel_err, imp );
    end
    
    if max(pres, dres) < params.opt_tol %|| imp < 1e-5
%     if pres < params.opt_tol
        break;
    end
    
    % update Lagrange multipliers
    for i = 1:length(V)
        V{i} = V{i} - tdiff{i}/mu;
    end
    
    if rem(iter, 10) == 0
        mu = max( mu*0.9, params.mu_min );
    end
end

results.X = Y;
results.E = E;
results.V = V;
results.T = data.T;
results.U = U;
results.iter = iter;
results.cpu = toc;
results.mu = mu;
% results.lambda = lambda;

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