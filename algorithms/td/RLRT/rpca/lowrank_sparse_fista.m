function results = lowrank_sparse_fista( data, params, funs )
% FISTA implementation for solving
%   min 0.5*||AY - T|| + \lambda_*\sum_i||Xi,(i)||_* + \lambda1*||E||_1
% where Z = [ X1,...,X_N,E ]', Y is the FISTA auxiliary variable, 
% A is a linear operator, T is given observation.

% tic;
% function handles
eval_f = funs.f;
eval_grad_f = funs.grad_f;
Aprod = funs.Aprod;
Atprod = funs.Atprod;
solve_prox_grad = funs.prox_grad_solver;

T = data.T;
if params.IsTC
    T = data.b;
end
N = length( size(params.X0) );
Y = cell(1,N+1);
for i = 1:N
    Y{i} = params.X0;
end
Y{N+1} = params.E0;
Z = Y;

use_cont = isfield( params, 'use_cont' ) && params.use_cont;
lambda1 = params.lambda;
lambdaS = params.lambdaS;
if use_cont
    if lambdaS > 0
        lambdaS_min = lambdaS;
        r = lambda1 / lambdaS;      lambdaS = 0.99*norm(T);     lambda1 = lambdaS*r;
    else
        r = params.rRatio /sqrt( max(size(data.T)) );
        lambdaS = 0.99*norm(T);     lambda1 = lambdaS*r;
        lambdaS_min = lambdaS * 1e-5;
    end
    eta = 0.97;
end

mu = params.mu0;
beta = 0.5;
t = 0;

fZ = eval_f( Z, Aprod, T );
trnorm = 0;
oneNorm = 0;

for iter = 1:params.max_iter
    
    % compute gradient at Y
    grad_f = eval_grad_f( Y, Aprod, Atprod, T );
    
    Zp = Z;     trnorm_p = trnorm;    fZp = fZ;     oneNorm_p = oneNorm;
    [ Z, trnorm ] = solve_prox_grad( Y, mu, grad_f, lambdaS, lambda1 );
    % backtrack linesearch
    fY = eval_f( Y, Aprod, T );
    ZYdiff = tensor_array_diff( Z, Y );
    approx = fY + tensor_array_innerprod( grad_f, ZYdiff ) + tensor_array_norm(ZYdiff)^2 / (2*mu);
    fZ = eval_f( Z, Aprod, T );
    while fZ > approx
        mu = mu * beta;
        [ Z, trnorm ] = solve_prox_grad( Y, mu, grad_f, lambdaS, lambda1 );
        ZYdiff = tensor_array_diff( Z, Y );
        approx = fY + tensor_array_innerprod( grad_f, ZYdiff ) + tensor_array_norm(ZYdiff)^2 / (2*mu);
        fZ = eval_f( Z, Aprod, T );
    end
    oneNorm = tensor_1norm(Z{N+1});
    
    % FISTA acceleration step
    tp = t;
    t = (1 + sqrt(1+4*tp^2))/2;
    Zdiff = tensor_array_diff( Z, Zp );
    Y = tensor_array_add( Zp, tensor_array_scale(Zdiff,(tp-1)/t) );
    
    
    % compute optimality stats
    Fp = fZp + lambdaS*trnorm_p + lambda1*oneNorm_p;
    F = fZ + lambdaS*trnorm + lambda1*oneNorm;
    rel_F = abs(F-Fp)/Fp;
    denom = tensor_array_norm(Zp);
    rel_X = tensor_array_norm(Zdiff) / denom;   if denom == 0; rel_X = 1; end
    rel_err = norm(ten_sum_all(Z(1:N))-data.X) / norm(data.X);
    
    % print
    if params.verbose && rem( iter, 50 ) == 0
        fprintf('Iter: %d,   fit: %3.2e,   rel_F: %3.2e,    rel_X: %3.2e,   rel_err: %3.2e, S: %3.2e\n', ...
            iter, fZ, rel_F, rel_X, rel_err, lambdaS);
    end
    
    if max(rel_F,rel_X) < params.opt_tol
        break;
    end
    
    % update lambda's
    if use_cont
        lambdaS = max( lambdaS*eta, lambdaS_min );
        lambda1 = lambdaS*r;
    end
end

results.vars = Z;
results.V = params.V0;
results.T = data.T;
results.iter = iter;
% results.cpu = toc;
results.mu = mu;
results.lambda1 = lambda1;
results.lambdaS = lambdaS;
results.IsTC = params.IsTC;

end