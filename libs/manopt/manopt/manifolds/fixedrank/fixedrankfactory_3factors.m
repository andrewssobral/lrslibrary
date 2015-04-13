function M = fixedrankfactory_3factors(m, n, k)
% Manifold of m-by-n matrices of rank k with polar quotient geometry.
%
% function M = fixedrankfactory_3factors(m, n, k)
%
% Follows the polar quotient geometry described in the following paper:
% G. Meyer, S. Bonnabel and R. Sepulchre,
% "Linear regression under fixed-rank constraints: a Riemannian approach",
% ICML 2011.
%
% Paper link: http://www.icml-2011.org/papers/350_icmlpaper.pdf
%
% Additional reference is
%
% B. Mishra, R. Meyer, S. Bonnabel and R. Sepulchre
% "Fixed-rank matrix factorizations and Riemannian low-rank optimization",
% arXiv, 2012.
%
% Paper link: http://arxiv.org/abs/1209.0430
%
% A point X on the manifold is represented as a structure with three
% fields: L, S and R. The matrices L (mxk) and R (nxk) are orthonormal,
% while the matrix S (kxk) is a symmetric positive definite full rank
% matrix.
%
% Tangent vectors are represented as a structure with three fields: L, S
% and R.

% This file is part of Manopt: www.manopt.org.
% Original author: Bamdev Mishra, Dec. 30, 2012.
% Contributors:
% Change log:
    
    M.name = @() sprintf('LSR'' quotient manifold of %dx%d matrices of rank %d', m, n, k);
    
    M.dim = @() (m+n-k)*k;
    
    % Choice of the metric on the orthnormal space is motivated by the symmetry present in the
    % space. The metric on the positive definite space is its natural metric.
    M.inner = @(X, eta, zeta) eta.L(:).'*zeta.L(:) + eta.R(:).'*zeta.R(:) ...
        + trace( (X.S\eta.S) * (X.S\zeta.S) );
    
    M.norm = @(X, eta) sqrt(M.inner(X, eta, eta));
    
    M.dist = @(x, y) error('fixedrankfactory_3factors.dist not implemented yet.');
    
    M.typicaldist = @() 10*k;
    
    skew = @(X) .5*(X-X');
    symm = @(X) .5*(X+X');
    stiefel_proj = @(L, H) H - L*symm(L'*H);
    
    M.egrad2rgrad = @egrad2rgrad;
    function eta = egrad2rgrad(X, eta)
        eta.L = stiefel_proj(X.L, eta.L);
        eta.S = X.S*symm(eta.S)*X.S;
        eta.R = stiefel_proj(X.R, eta.R);
    end
    
    
    M.ehess2rhess = @ehess2rhess;
    function Hess = ehess2rhess(X, egrad, ehess, eta)
        
        % Riemannian gradient for the factor S
        rgrad.S = X.S*symm(egrad.S)*X.S;
        
        % Directional derivatives of the Riemannian gradient
        Hess.L = ehess.L - eta.L*symm(X.L'*egrad.L);
        Hess.L = stiefel_proj(X.L, Hess.L);
        
        Hess.R = ehess.R - eta.R*symm(X.R'*egrad.R);
        Hess.R = stiefel_proj(X.R, Hess.R);
        
        Hess.S = X.S*symm(ehess.S)*X.S +  2*symm(eta.S*symm(egrad.S)*X.S);
        
        % Correction factor for the non-constant metric on the factor S
        Hess.S = Hess.S - symm(eta.S*(X.S\rgrad.S));
        
        % Projection onto the horizontal space
        Hess = M.proj(X, Hess);
    end
    
    
    M.proj = @projection;
    function etaproj = projection(X, eta)
        % First, projection onto the tangent space of the total sapce
        eta.L = stiefel_proj(X.L, eta.L);
        eta.R = stiefel_proj(X.R, eta.R);
        eta.S = symm(eta.S);
        
        % Then, projection onto the horizontal space
        SS = X.S*X.S;
        AS = X.S*(skew(X.L'*eta.L) + skew(X.R'*eta.R) - 2*skew(X.S\eta.S))*X.S;
        omega = lyap(SS, -AS);
        
        etaproj.L = eta.L - X.L*omega;
        etaproj.S = eta.S - (X.S*omega - omega*X.S);
        etaproj.R = eta.R - X.R*omega;
    end
    
    M.tangent = M.proj;
    M.tangent2ambient = @(X, eta) eta;
    
    M.retr = @retraction;
    function Y = retraction(X, eta, t)
        if nargin < 3
            t = 1.0;
        end
        
        L = chol(X.S);
        Y.S = L'*expm(L'\(t*eta.S)/L)*L;
        Y.L = uf(X.L + t*eta.L);
        Y.R = uf(X.R + t*eta.R);
    end
    
    M.exp = @exponential;
    function Y = exponential(X, eta, t)
        if nargin < 3
            t = 1.0;
        end
        Y = retraction(X, eta, t);
        warning('manopt:fixedrankfactory_3factors:exp', ...
            ['Exponential for fixed rank ' ...
            'manifold not implemented yet. Lsed retraction instead.']);
    end
    
    M.hash = @(X) ['z' hashmd5([X.L(:) ; X.S(:) ; X.R(:)])];
    
    M.rand = @random;
    % Factors L and R live on Stiefel manifolds, hence we will reuse
    % their random generator.
    stiefelm = stiefelfactory(m, k);
    stiefeln = stiefelfactory(n, k);
    function X = random()
        X.L = stiefelm.rand();
        X.R = stiefeln.rand();
        X.S = diag(1+rand(k, 1));
    end
    
    M.randvec = @randomvec;
    function eta = randomvec(X)
        % A random vector on the horizontal space
        eta.L = randn(m, k);
        eta.R = randn(n, k);
        eta.S = randn(k, k);
        eta = projection(X, eta);
        nrm = M.norm(X, eta);
        eta.L = eta.L / nrm;
        eta.R = eta.R / nrm;
        eta.S = eta.S / nrm;
    end
    
    M.lincomb = @lincomb;
    
    M.zerovec = @(X) struct('L', zeros(m, k), 'S', zeros(k, k), ...
        'R', zeros(n, k));
    
    M.transp = @(x1, x2, d) projection(x2, d);
    
    % vec and mat are not isometries, because of the unusual inner metric.
    M.vec = @(X, U) [U.L(:) ; U.S(:); U.R(:)];
    M.mat = @(X, u) struct('L', reshape(u(1:(m*k)), m, k), ...
        'S', reshape(u((m*k+1): m*k + k*k), k, k), ...
        'R', reshape(u((m*k+ k*k + 1):end), n, k));
    M.vecmatareisometries = @() false;
    
end

% Linear combination of tangent vectors
function d = lincomb(x, a1, d1, a2, d2) %#ok<INLSL>
    
    if nargin == 3
        d.L = a1*d1.L;
        d.R = a1*d1.R;
        d.S = a1*d1.S;
    elseif nargin == 5
        d.L = a1*d1.L + a2*d2.L;
        d.R = a1*d1.R + a2*d2.R;
        d.S = a1*d1.S + a2*d2.S;
    else
        error('Bad use of fixedrankfactory_3factors.lincomb.');
    end
    
end

function A = uf(A)
    [L, unused, R] = svd(A, 0); %#ok
    A = L*R';
end
