function [  x, w ] =  solveWXADMM( U, v, param)


    if isfield(param,'TOL'),
        tol = param.tol;
    else
        tol = 1e-4;
    end

    if isfield(param,'maxIter'),
        maxIter = param.maxIter;
    else
        maxIter = 20;
    end

    %% Data preprocessing

    [m, n] = size(U);
    x = zeros(m,1);


    %% ADMM solver

    G = param.admm.G; % G is  a n*g matrix, with each column correspobding to one group

    Z = param.admm.Z;

    Y = param.admm.Y; % Y is 

    lambda =  param.lambda;
    mu =  0.3/mean(abs(v));

    gamma = 1.6;

    x0 = x;
    Im = eye(param.rank);

    for iter=1: maxIter,


        A = [lambda * Im, lambda * U'; lambda *U,  lambda *  param.admm.I +  mu*  diag(sum(G,2))  ];

        b = [lambda *U'*v; lambda*v - sum(G.*Y, 2) + mu * sum(G.*Z, 2) ];

        b=sparse(b);
        
        wx = A \ b;
  
        x = wx(param.rank+1:end);

        % Z update    
        X=diag(sparse(x));

        DX = X*G;  % DX: each column is grouped varianbles

        R = DX + Y/mu;
 
        Z = shrinkageMex(R,   [1, mu] );

        %   Y update 
        Y = Y +    mu * (   DX - Z);   

        mu = mu*gamma;

        xDiff = norm(x-x0);

        if (xDiff < tol),
    %         disp('**********************Converged*******************');
            break;   
        end
        x0 = x;
    end

    w = wx(1:param.rank);
    x=full(x);


end

