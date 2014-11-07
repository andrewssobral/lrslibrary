function results = rpca_for_tensor( data, params )
% params.mode = 0: apply RPCA to each frontal slice of the tensor. only
% supported for IsTC = false.
% params.mode = i: apply RPCA to the i-th mode unfolding of the tensor
% Uses IALM by Yi Ma's group
% For partial observations case, RPCA does not support mode 0 yet.

tic;
if ~params.IsTC
switch params.mode
    case 0
        iters = [];
        for k = 1:size(data.T,3)
            T = double( data.T(:,:,k) );
            lambda = 1 / sqrt(max(size(T)));
            [Xk Ek iter] = inexact_alm_rpca(T, lambda*params.rRatio, params.opt_tol, params.max_iter, 1/(params.mu1fac*std(T(:))) );
            results.X(:,:,k) = Xk;
            results.E(:,:,k) = Ek;
            iters(k) = iter;
        end
        results.iter = mean(iters);
    otherwise
%     N = length(size(double(data.T)));
        Tmat = tenmat(data.T,params.mode);
        T = double( Tmat );
        lambda = 1 / sqrt(max(size(T)));
        [X E iter] = inexact_alm_rpca(T, lambda*params.rRatio, params.opt_tol, params.max_iter, 1/params.mu1 );
        results.X = tensor(tenmat( X, Tmat.rdims, Tmat.cdims, Tmat.tsize) );
        results.E = tensor(tenmat( E, Tmat.rdims, Tmat.cdims, Tmat.tsize) );
        results.iter = iter;
end

else
    % RPCA with missing data
    results = rpca_tc( data, params );
end

results.cpu = toc;
results.IsTC = params.IsTC;

end

% ======================================================================= %
function results = rpca_tc( data, params )
% solves
%   min_{X,E} ||X||_* + \lambda_1||E||_1
%   s.t.      Ac(Y+E) = Tc
%             Y = X
%
% data.matInd: indices in Omega w.r.t. i-th mode of T

Tmat = tenmat(data.T,params.mode);
Xmat = double(tenmat(data.X,params.mode));
[ n, m ] = size(double(Tmat));

data = tenInd2matInd_core( data, params.mode );
Omega = data.matInd;
b = double(Tmat(Omega));
X = zeros( n, m );
Y = X;
E = zeros( n, m );

%%%%%%%% Fix lambda value here!! %%%%%%%%%%
% lambda = params.lambda;
lambda = params.rRatio / sqrt(max(m,n));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = params.mu1;

V = cell(1,2);
V{1} = zeros( size(b) );
V{2} = zeros( size(X) );
bnorm2 = norm(b)^2;

% assign operators
Aprod = @(X_param)A_select( X_param, Omega );
Atprod = @(b_param)At_select( b_param, Omega, size(E) );

for iter = 1:params.max_iter
    % solve X_i's
    P = Y - mu*V{2};
    X = matrix_shrinkage( P, mu );
    
    % solve E
    p = b + mu*V{1} - Aprod(Y);     P = Atprod(p);
    E(:) = shrinkage_v( P(:), mu*lambda );
    
    % solve Y
    Yp = Y;
    Y = X + mu*V{2};
    Y(Omega) = (Y(Omega) - E(Omega) + b + mu*V{1}) / 2;
    
    % compute optimality stats
    tdiff{1} = Aprod( Y + E ) - b;
    tdiff{2} = Y - X;
    Ydiff = Y - Yp;
    YdiffOmega = Ydiff(Omega);
    YpOmega = Yp(Omega);
    
    pres = norm( [tdiff{1}; tdiff{2}(:)] ) / norm( [b; X(:)] );
    dres = norm( [YdiffOmega(:); Ydiff(:)] ) / norm( [YpOmega(:); Y(:)] );
    
    rel_err = norm(Y-Xmat, 'fro') / norm(Xmat,'fro');
    
    % print
    if params.verbose
        fprintf('Iter: %d,   pinf: %3.2e,   dinf: %3.2e,    rel_err: %3.2e\n', iter, pres, dres, rel_err );
    end
    
    if max(pres, dres) < params.opt_tol
%     if pres < params.opt_tol
        break;
    end
    
    % update Lagrange multipliers
    V{1} = V{1} - tdiff{1}/mu;
    V{2} = V{2} - tdiff{2}/mu;
    
%     mu = max(mu / 1.5, params.mu_min);
end

results.X = tensor( tenmat( X, Tmat.rdims, Tmat.cdims, Tmat.tsize ) );
results.E = tensor( tenmat( E, Tmat.rdims, Tmat.cdims, Tmat.tsize ) );
results.V = V;
results.b = b;
results.iter = iter;
results.mu = mu;
results.lambda = lambda;


end


function AX = A_select( X, Omega )
% AX is a vector

AX = X(Omega);

end

function X = At_select( b, Omega, sz )

X = zeros( sz );
X(Omega) = b;
end