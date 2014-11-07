function results = tensor_rpca_adal( data, params )
% Solve
%       min_{X,E} \sum_i ||X_(i)||_* + \lambda1*||E||_1
%       s.t. X + E = T

% data.T
% params
% X, V are cell arrays of tensors.
%
% Algorithm: ADAL

tic;
T = data.T;
N = length( size(params.X0) );
X = cell( 1, N );
U = cell( 1, N );
V = cell( 1, N+1 );
for i = 1:N
    X{i} = params.X0;
    V{i} = params.V0{i};
end
V{N+1} = params.V0{1};
Y = params.X0;
E = params.E0;

lambda = params.lambda;
mu1 = params.mu1;
mu2 = params.mu2;
verbose = params.verbose;


for iter = 1:params.max_iter
    % solve X_i's
    for i = 1:N
        %%% for PROPACK %%%
%         tmode = i;
        %%%%%%%%%%%%%%%%%%% 
        [X{i}, junk, U{i}] = tensor_shrinkage( Y+mu2*V{i}, mu2, i );
    end
    
    % solve E
    P = T - Y + mu1*V{N+1};
    E = shrinkage_t( P, lambda*mu1 );
    
    % solve Y
    Yrhs = V{N+1} + (T - E)/mu1;
    Yrhs = Yrhs + ten_sum_all( X )/mu2 - ten_sum_all( V(1:N) );
    Yprev = Y;
    Y = Yrhs / (N/mu2 + 1/mu1);
    
    % compute optimality stats
    pres = 0;
    tdiff = cell( 1, N+1 );
    for i = 1:N
        tdiff{i} = X{i} - Y;
        pres = pres + norm( tdiff{i} )^2;
    end
    tdiff{N+1} = Y + E - T;
    pres = pres + norm(tdiff{N+1})^2;     %pres = sqrt(pres);
    
    pres = sqrt( pres / (norm(T)^2+N*norm(Y)^2) );
    Ydiff = Y-Yprev;
    dres = norm(Ydiff) / norm(Yprev);
    rel_err = norm(Y-data.X) / norm(data.X);
    
    % print
    if verbose
    fprintf('Iter: %d,   pinf: %3.2e,   dinf: %3.2e,    rel_err: %3.2e\n', iter, pres, dres, rel_err );
    end
    
    if max(pres, dres) < params.opt_tol
%     if pres < params.opt_tol
        break;
    end
    
    % update Lagrange multipliers
    for i = 1:N
        V{i} = V{i} - tdiff{i}/mu2;
    end
    V{N+1} = V{N+1} - tdiff{N+1}/mu1;
    
%     if rem(iter,50) == 0
%         mu1 = max(mu1*0.8, params.mu_min);    mu2 = max(mu2*0.8, params.mu_min);
%     end
end

results.X = Y;
results.E = E;
results.V = V;
results.T = T;
results.U = U;
results.iter = iter;
results.cpu = toc;
results.mu = mu1;
results.lambda = lambda;

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