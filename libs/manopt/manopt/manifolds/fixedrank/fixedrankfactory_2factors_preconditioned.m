function M = fixedrankfactory_2factors_preconditioned(m, n, k)
% Manifold of m-by-n matrices of rank k with new balanced quotient geometry
%
% function M = fixedrankfactory_2factors_preconditioned(m, n, k)
%
% This follows the quotient geometry described in the following paper:
% B. Mishra, K. Adithya Apuroop and R. Sepulchre,
% "A Riemannian geometry for low-rank matrix completion",
% arXiv, 2012.
%
% Paper link: http://arxiv.org/abs/1211.1550
%
% This geoemtry is tuned to least square problems such as low-rank matrix
% completion.
%
% A point X on the manifold is represented as a structure with two
% fields: L and R. The matrices L (mxk) and R (nxk) are full column-rank
% matrices.
%
% Tangent vectors are represented as a structure with two fields: L, R

% This file is part of Manopt: www.manopt.org.
% Original author: Bamdev Mishra, Dec. 30, 2012.
% Contributors:
% Change log:
    
    
    
    M.name = @() sprintf('LR''(tuned for least square problems) quotient manifold of %dx%d matrices of rank %d', m, n, k);
    
    M.dim = @() (m+n-k)*k;
    
    
    % Some precomputations at the point X to be used in the inner product (and
    % pretty much everywhere else).
    function X = prepare(X)
        if ~all(isfield(X,{'LtL','RtR','invRtR','invLtL'}))
            L = X.L;
            R = X.R;
            X.LtL = L'*L;
            X.RtR = R'*R;
            X.invLtL = inv(X.LtL);
            X.invRtR = inv(X.RtR);
        end
    end
    
    
    % The choice of metric is motivated by symmetry and tuned to least square
    % objective function
    M.inner = @iproduct;
    function ip = iproduct(X, eta, zeta)
        X = prepare(X);
        
        ip = trace(X.RtR*(eta.L'*zeta.L)) + trace(X.LtL*(eta.R'*zeta.R));
    end
    
    M.norm = @(X, eta) sqrt(M.inner(X, eta, eta));
    
    M.dist = @(x, y) error('fixedrankfactory_2factors_preconditioned.dist not implemented yet.');
    
    M.typicaldist = @() 10*k;
    
    symm = @(M) .5*(M+M');
    
    M.egrad2rgrad = @egrad2rgrad;
    function eta = egrad2rgrad(X, eta)
        X = prepare(X);
        
        eta.L = eta.L*X.invRtR;
        eta.R = eta.R*X.invLtL;
    end
    
    M.ehess2rhess = @ehess2rhess;
    function Hess = ehess2rhess(X, egrad, ehess, eta)
        X = prepare(X);
        
        % Riemannian gradient
        rgrad = egrad2rgrad(X, egrad);
        
        % Directional derivative of the Riemannian gradient
        Hess.L = ehess.L*X.invRtR - 2*egrad.L*(X.invRtR * symm(eta.R'*X.R) * X.invRtR);
        Hess.R = ehess.R*X.invLtL - 2*egrad.R*(X.invLtL * symm(eta.L'*X.L) * X.invLtL);
        
        % We still need a correction factor for the non-constant metric
        Hess.L = Hess.L + rgrad.L*(symm(eta.R'*X.R)*X.invRtR) + eta.L*(symm(rgrad.R'*X.R)*X.invRtR) - X.L*(symm(eta.R'*rgrad.R)*X.invRtR);
        Hess.R = Hess.R + rgrad.R*(symm(eta.L'*X.L)*X.invLtL) + eta.R*(symm(rgrad.L'*X.L)*X.invLtL) - X.R*(symm(eta.L'*rgrad.L)*X.invLtL);
        
        % Project on the horizontal space
        Hess = M.proj(X, Hess);
        
    end
    
    M.proj = @projection;
    function etaproj = projection(X, eta)
        X = prepare(X);
        
        Lambda =  (eta.R'*X.R)*X.invRtR  -   X.invLtL*(X.L'*eta.L);
        Lambda = Lambda/2;
        etaproj.L = eta.L + X.L*Lambda;
        etaproj.R = eta.R - X.R*Lambda';
    end
    
    M.tangent = M.proj;
    M.tangent2ambient = @(X, eta) eta;
    
    
    
    M.retr = @retraction;
    function Y = retraction(X, eta, t)
        if nargin < 3
            t = 1.0;
        end
        Y.L = X.L + t*eta.L;
        Y.R = X.R + t*eta.R;
        
        % Numerical conditioning step: A simpler version.
        % We need to ensure that L and R are do not have very relative
        % skewed norms.
        
        scaling = norm(X.L, 'fro')/norm(X.R, 'fro');
        scaling = sqrt(scaling);
        Y.L = Y.L / scaling;
        Y.R = Y.R * scaling;
        
        % These are reused in the computation of the gradient and Hessian
        Y = prepare(Y);
    end
    
    
    M.exp = @exponential;
    function Y = exponential(X, eta, t)
        if nargin < 3
            t = 1.0;
        end
        
        Y = retraction(X, eta, t);
        warning('manopt:fixedrankfactory_2factors_preconditioned:exp', ...
            ['Exponential for fixed rank ' ...
            'manifold not implemented yet. Used retraction instead.']);
    end
    
    M.hash = @(X) ['z' hashmd5([X.L(:) ; X.R(:)])];
    
    M.rand = @random;
    
    function X = random()
        X.L = randn(m, k);
        X.R = randn(n, k);
    end
    
    M.randvec = @randomvec;
    function eta = randomvec(X)
        eta.L = randn(m, k);
        eta.R = randn(n, k);
        eta = projection(X, eta);
        nrm = M.norm(X, eta);
        eta.L = eta.L / nrm;
        eta.R = eta.R / nrm;
    end
    
    M.lincomb = @lincomb;
    
    M.zerovec = @(X) struct('L', zeros(m, k),'R', zeros(n, k));
    
    M.transp = @(x1, x2, d) projection(x2, d);
    
    % vec and mat are not isometries, because of the unusual inner metric.
    M.vec = @(X, U) [U.L(:) ; U.R(:)];
    M.mat = @(X, u) struct('L', reshape(u(1:(m*k)), m, k), ...
        'R', reshape(u((m*k+1):end), n, k));
    M.vecmatareisometries = @() false;
    
end

% Linear combination of tangent vectors
function d = lincomb(x, a1, d1, a2, d2) %#ok<INUSL>
    
    if nargin == 3
        d.L = a1*d1.L;
        d.R = a1*d1.R;
    elseif nargin == 5
        d.L = a1*d1.L + a2*d2.L;
        d.R = a1*d1.R + a2*d2.R;
    else
        error('Bad use of fixedrankfactory_2factors_preconditioned.lincomb.');
    end
    
end





