function M = fixedrankfactory_2factors(m, n, k)
% Manifold of m-by-n matrices of rank k with balanced quotient geometry.
%
% function M = fixedrankfactory_2factors(m, n, k)
%
% This follows the balanced quotient geometry described in the following paper:
% G. Meyer, S. Bonnabel and R. Sepulchre,
% "Linear regression under fixed-rank constraints: a Riemannian approach",
% ICML 2011.
%
% Paper link: http://www.icml-2011.org/papers/350_icmlpaper.pdf
%
% A point X on the manifold is represented as a structure with two
% fields: L and R. The matrices L (mxk) and R (nxk) are full column-rank
% matrices such that X = L*R'.
%
% Tangent vectors are represented as a structure with two fields: L, R

% This file is part of Manopt: www.manopt.org.
% Original author: Bamdev Mishra, Dec. 30, 2012.
% Contributors:
% Change log:
%   July 10, 2013 (NB) : added vec, mat, tangent, tangent2ambient
    
    
    M.name = @() sprintf('LR'' quotient manifold of %dx%d matrices of rank %d', m, n, k);
    
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
    
    % Choice of the metric is motivated by the symmetry present in the space
    M.inner = @iproduct;
    function ip = iproduct(X, eta, zeta)
        X = prepare(X);
        ip = trace(X.invLtL*(eta.L'*zeta.L)) + trace( X.invRtR*(eta.R'*zeta.R));
    end
    
    M.norm = @(X, eta) sqrt(M.inner(X, eta, eta));
    
    M.dist = @(x, y) error('fixedrankfactory_2factors.dist not implemented yet.');
    
    M.typicaldist = @() 10*k;
    
    symm = @(M) .5*(M+M');
    
    M.egrad2rgrad = @egrad2rgrad;
    function eta = egrad2rgrad(X, eta)
        X = prepare(X);
        eta.L = eta.L*X.LtL;
        eta.R = eta.R*X.RtR;
    end
    
    M.ehess2rhess = @ehess2rhess;
    function Hess = ehess2rhess(X, egrad, ehess, eta)
        X = prepare(X);
        
        % Riemannian gradient
        rgrad = egrad2rgrad(X, egrad);
        
        % Directional derivative of the Riemannian gradient
        Hess.L = ehess.L*X.LtL + 2*egrad.L*symm(eta.L'*X.L);
        Hess.R = ehess.R*X.RtR + 2*egrad.R*symm(eta.R'*X.R);
        
        % We need a correction term for the non-constant metric
        Hess.L = Hess.L - rgrad.L*((X.invLtL)*symm(X.L'*eta.L)) - eta.L*(X.invLtL*symm(X.L'*rgrad.L)) + X.L*(X.invLtL*symm(eta.L'*rgrad.L));
        Hess.R = Hess.R - rgrad.R*((X.invRtR)*symm(X.R'*eta.R)) - eta.R*(X.invRtR*symm(X.R'*rgrad.R)) + X.R*(X.invRtR*symm(eta.R'*rgrad.R));
        
        % Projection onto the horizontal space
        Hess = M.proj(X, Hess);
    end
    
    M.proj = @projection;
    % Projection of the vector eta onto the horizontal space
    function etaproj = projection(X, eta)
        X = prepare(X);
        
        SS = (X.LtL)*(X.RtR);
        AS = (X.LtL)*(X.R'*eta.R) - (eta.L'*X.L)*(X.RtR);
        Omega = lyap(SS, SS,-AS);
        etaproj.L = eta.L + X.L*Omega';
        etaproj.R = eta.R - X.R*Omega;
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
        % We need to ensure that L and R do not have very relative
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
        warning('manopt:fixedrankfactory_2factors:exp', ...
            ['Exponential for fixed rank ' ...
            'manifold not implemented yet. Used retraction instead.']);
    end
    
    M.hash = @(X) ['z' hashmd5([X.L(:) ; X.R(:)])];
    
    M.rand = @random;
    function X = random()
        % A random point on the total space
        X.L = randn(m, k);
        X.R = randn(n, k);
        X = prepare(X);
    end
    
    M.randvec = @randomvec;
    function eta = randomvec(X)
        % A random vector in the horizontal space
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
        error('Bad use of fixedrankfactory_2factors.lincomb.');
    end
    
end





