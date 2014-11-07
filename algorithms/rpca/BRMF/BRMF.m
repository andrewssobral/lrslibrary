function [U V Tau LowRank] = BRMF(Y, opts, U, V)
    rand('state',0);
    randn('state',0);
    
    
    [m n] = size(Y);
    if nargin < 3
       U = randn(m, opts.r);
       V = randn(n, opts.r);
    end
    
    accU = zeros(m, opts.r);
    accV = zeros(n, opts.r);
    accTau = zeros(m, n);
    LowRank = zeros(m, n);
    lambda = ones(m, n);
    for iter = 1 : opts.maxIter
         disp(iter);
        %% Sample hyperparamters
        [mu_u Lambda_u] = sampleHyper(U, m, opts.invW_0, opts.beta_0, opts.nu_0, opts.r);
        [mu_v Lambda_v] = sampleHyper(V, n, opts.invW_0, opts.beta_0, opts.nu_0, opts.r);
        %% Sample Tau
        res = abs(Y - U * V') ;
        Tau = drawFromIG(sqrt(lambda) ./ res, lambda);
        Tau = 1 ./ Tau;
        Tau = Tau + 1e-6;
        
       %% Sample lambda
        
        lambda = drawFromIG(sqrt((Tau + opts.a) ./ opts.b), Tau + opts.a);
        lambda = 1 ./ lambda;
        %% Sample U
        parfor i = 1 : m
            Lambda_i = V' * spdiags(1 ./ Tau(i, :)', 0, n, n) * V;
            Lambda_i = inv(Lambda_u + Lambda_i);
            u_i = Lambda_i * (sum(bsxfun(@times, (Y(i, :) ./ Tau(i, :))', V))'  + Lambda_u * mu_u);
            lam = chol(Lambda_i); lam=lam'; 
            U(i, :) = lam * randn(opts.r, 1) + u_i;
            
        end
        
        %% Sample V
        parfor i = 1 : n
            Lambda_i = U' * spdiags(1 ./ Tau(:, i), 0, m, m) * U;
            Lambda_i = inv(Lambda_v + Lambda_i);
            Lambda_i = (Lambda_i + Lambda_i') / 2;
            v_i = Lambda_i * (sum(bsxfun(@times, Y(:, i) ./ Tau(:, i), U))' + Lambda_v * mu_v);
            lam = chol(Lambda_i); lam=lam'; 
            V(i, :) = lam * randn(opts.r, 1) + v_i;
        end
        
        if (iter > opts.burnin)
            accU = accU + U;
            accV = accV + V;
            accTau = accTau + Tau;
            LowRank = LowRank + U * V';
        end
        
    end
    U = accU / (opts.maxIter - opts.burnin);
    V = accV / (opts.maxIter - opts.burnin);
    Tau = accTau / (opts.maxIter - opts.burnin);
    LowRank = LowRank / (opts.maxIter - opts.burnin);
end