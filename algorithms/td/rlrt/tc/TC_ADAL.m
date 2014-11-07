function results = TC_ADAL( data, params )

tic;
Y = data.X;
b = data.b;
Omega = data.Omega;

mu = params.mu0;
N = length( size(Y) );
lambda = params.lambda;

Lamb = cell( 1, N );
Xs = cell(1,N);
for i = 1:N
    Lamb{i} = data.Lamb{i};
    Xs{i} = data.X;
end

for iter = 1:params.max_iter
    % solve for X_{N+1} = Y
    Yprev = Y;
    R = tenzeros( size(Y) );
    R(Omega) = lambda*b;
    R = R + ten_sum_all( Lamb ) + ten_sum_all( Xs )/mu;
    Y = R * mu/N;
    Y(Omega) = R(Omega) / (lambda+N/mu);
    
    % solve for X_i's
%     Xs = cell(1,N);
    for i = 1:N
        Xs{i} = tensor_shrinkage( Y-mu*Lamb{i}, mu, i );
        Lamb{i} = Lamb{i} - (Y-Xs{i})/mu;
    end
    
    % compute optimality stats
    pinf = 0;
    for i = 1:N
        pinf = pinf + norm( double(tenmat( Xs{i}-Y, i )), 'fro' )^2;
    end
    ynorm = norm(double(tenmat(Y,1)),'fro');
    pinf = sqrt(pinf / N) / ynorm;
%     pinf = pinf / ynorm;
    dinf = norm( double(tenmat(Y - Yprev, 1 )), 'fro' ) / ynorm;
%     dinf = dinf / ynorm;
    
    % print
    fprintf('Iter: %d,   pinf: %3.2e,   dinf: %3.2e\n', iter, pinf, dinf);
    
    if pinf < params.opt_tol && dinf < params.opt_tol
        break;
    end
    
    % update mu
%     if rem( iter, 10 ) == 0
%         mu = mu / 5;
%     end
end

diffNorm = norm( double( tenmat( Y - data.X, 1 )), 'fro' );
trueNorm = norm( double( tenmat( data.X, 1 )), 'fro' );
rel_err = diffNorm / trueNorm;

% save results
results.X = Y;
results.iter = iter;
results.cpu = toc;
results.rel_err = rel_err;

end