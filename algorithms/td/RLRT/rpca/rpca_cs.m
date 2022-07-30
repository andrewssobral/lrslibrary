function results = rpca_cs( data, params )
% Solve
%       min_{X,E} \sum_i ||X_(i)||_* + \lambda1*||E||_1
%       s.t. A_Q(F'*X + E) = T_Q
% A_Q = partial Fourier transform, F = DFT
% Reformulate the problem as
% 
%       min_{Xi's, E} \sum_i ||Xi_(i)||_* + \lambda1*||E||_1
%       s.t. A_Omega(Xi) + A_Q(E) = T_Q, i = 1,...,N.
%
% data.T
% params
% X, V are cell arrays of tensors.
%
% Algorithm: I-ADAL

tic;
b = data.b;
% N = ndims(params.X0);
X = params.X0;
% U = cell( 1, N );
V = zeros( size(b) );
E = params.E0;

lambda1 = params.lambda;
mu = 10;    %params.mu1;
eta = 1/1.1;    %params.eta;
verbose = params.verbose;
bnorm = norm(b);

% function handles
Aprod = data.funs.Aprod;
Atprod = data.funs.Atprod;
Asel = data.funs.Aprod;
Atsel = data.funs.Atprod;
% Asel = data.funs.Asel;
% Atsel = data.funs.Atsel;


for iter = 1:params.max_iter
    % solve X_i's
    gradX = Atsel( Asel(X) + Aprod(E) - b - mu*V );
    PX = X - eta*gradX;
    [X, junk, U] = tensor_shrinkage( PX, mu*eta, 4 );
    
    % solve E
    Ep = E;
    D = Asel(X) - b - mu*V;
    gradE = Atprod( Aprod(E) + D );
    PE = E - eta*gradE;
    E = shrinkage_t( PE, lambda1*mu*eta );
    
    % compute optimality stats
    pres = 0;
    tdiff = Asel(X) + Aprod(E) - b;
    pres = norm(tdiff);
    
    Ediff = E - Ep;
    pres = pres / bnorm;
    dres = norm( Atsel(Aprod(Ediff)) ) / norm( Atsel(Aprod(Ep)) );
    rel_err = norm(X-data.X) / norm(data.X);
    
    % print
    if verbose
    fprintf('Iter: %d,   pinf: %3.2e,   dinf: %3.2e,    rel_err: %3.2e\n', iter, pres, dres, rel_err );
    end
    
    if max(pres, dres) < params.opt_tol
        break;
    end
    
    % update Lagrange multipliers
    V = V - tdiff/mu;
    
%     if rem(iter,10) == 0
%         mu = max(mu*0.9, params.mu_min);
%     end
end

results.X = X;  %real(iffn(Y));
results.E = E;
results.V = V;
results.T = b;
results.U = U;
results.iter = iter;
results.cpu = toc;
results.mu = mu;
results.lambda1 = lambda1;
results.eta = eta;
end
