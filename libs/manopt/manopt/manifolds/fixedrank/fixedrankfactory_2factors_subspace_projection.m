function M = fixedrankfactory_2factors_subspace_projection(m, n, k)
% Manifold of m-by-n matrices of rank k with quotient geometry.
%
% function M = fixedrankfactory_2factors_subspace_projection(m, n, k)
%
% This follows the quotient geometry described in the following paper:
% B. Mishra, G. Meyer, S. Bonnabel and R. Sepulchre
% "Fixed-rank matrix factorizations and Riemannian low-rank optimization",
% arXiv, 2012.
%
% Paper link: http://arxiv.org/abs/1209.0430
%
% A point X on the manifold is represented as a structure with two
% fields: L and R. The matrices L (mxk) is orthonormal,
% while the matrix R (nxk) is a full column-rank
% matrix.
%
% Tangent vectors are represented as a structure with two fields: L, R.

% This file is part of Manopt: www.manopt.org.
% Original author: Bamdev Mishra, Dec. 30, 2012.
% Contributors:
% Change log:
    
    
    
    M.name = @() sprintf('LR'' quotient manifold of %dx%d matrices of rank %d', m, n, k);
    
    M.dim = @() (m+n-k)*k;
    
    % Some precomputations at the point X to be used in the inner product (and
    % pretty much everywhere else).
    function X = prepare(X)
        if ~all(isfield(X,{'RtR','invRtR'}) == 1)
            X.RtR = X.R'*X.R;
            X.invRtR = eye(size(X.R,2))/ X.RtR;
        end
    end
    
    % The choice of the metric is motivated by symmetry and scale
    % invariance in the total space
    M.inner = @iproduct;
    function ip = iproduct(X, eta, zeta)
        X = prepare(X);
        
        ip = eta.L(:).'*zeta.L(:)  + trace(X.invRtR*(eta.R'*zeta.R) );
    end
    
    M.norm = @(X, eta) sqrt(M.inner(X, eta, eta));
    
    M.dist = @(x, y) error('fixedrankfactory_2factors_subspace_projection.dist not implemented yet.');
    
    M.typicaldist = @() 10*k;
    
    skew = @(X) .5*(X-X');
    symm = @(X) .5*(X+X');
    stiefel_proj = @(L, H) H - L*symm(L'*H);
    
    M.egrad2rgrad = @egrad2rgrad;
    function eta = egrad2rgrad(X, eta)
        X = prepare(X);
        
        eta.L = stiefel_proj(X.L, eta.L);
        eta.R = eta.R*X.RtR;
    end
    
    
    M.ehess2rhess = @ehess2rhess;
    function Hess = ehess2rhess(X, egrad, ehess, eta)
        X = prepare(X);
        
        % Riemannian gradient
        rgrad = egrad2rgrad(X, egrad);
        
        % Directional derivative of the Riemannian gradient
        Hess.L = ehess.L - eta.L*symm(X.L'*egrad.L);
        Hess.L = stiefel_proj(X.L, Hess.L);
        
        Hess.R = ehess.R*X.RtR + 2*egrad.R*symm(eta.R'*X.R);
        
        % Correction factor for the non-constant metric on the factor R
        Hess.R = Hess.R - rgrad.R*((X.invRtR)*symm(X.R'*eta.R)) - eta.R*(X.invRtR*symm(X.R'*rgrad.R)) + X.R*(X.invRtR*symm(eta.R'*rgrad.R));
        
        % Projection onto the horizontal space
        Hess = M.proj(X, Hess);
    end
    
    
    M.proj = @projection;
    function etaproj = projection(X, eta)
        X = prepare(X);
        
        eta.L = stiefel_proj(X.L, eta.L); % On the tangent space
        SS = X.RtR;
        AS1 = 2*X.RtR*skew(X.L'*eta.L)*X.RtR;
        AS2 = 2*skew(X.RtR*(X.R'*eta.R));
        AS  = skew(AS1 + AS2);
        
        Omega = nested_sylvester(SS,AS);
        etaproj.L = eta.L - X.L*Omega;
        etaproj.R = eta.R - X.R*Omega;
    end
    
    M.tangent = M.proj;
    M.tangent2ambient = @(X, eta) eta;
    
    M.retr = @retraction;
    function Y = retraction(X, eta, t)
        if nargin < 3
            t = 1.0;
        end
        Y.L = uf(X.L + t*eta.L);
        Y.R = X.R + t*eta.R;
        
        % These are reused in the computation of the gradient and Hessian
        Y = prepare(Y);
    end
    
    M.exp = @exponential;
    function R = exponential(X, eta, t)
        if nargin < 3
            t = 1.0;
        end
        
        R = retraction(X, eta, t);
        warning('manopt:fixedrankfactory_2factors_subspace_projection:exp', ...
            ['Exponential for fixed rank ' ...
            'manifold not implemented yet. Lsed retraction instead.']);
    end
    
    M.hash = @(X) ['z' hashmd5([X.L(:) ; X.R(:)])];
    
    M.rand = @random;
    % Factors L lives on Stiefel manifold, hence we will reuse
    % its random generator.
    stiefelm = stiefelfactory(m, k);
    function X = random()
        X.L = stiefelm.rand();
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
    
    M.zerovec = @(X) struct('L', zeros(m, k),...
        'R', zeros(n, k));
    
    M.transp = @(x1, x2, d) projection(x2, d);
    
    % vec and mat are not isometries, because of the unusual inner metric.
    M.vec = @(X, U) [U.L(:) ; U.R(:)];
    M.mat = @(X, u) struct('L', reshape(u(1:(m*k)), m, k), ...
        'R', reshape(u((m*k+1):end), n, k));
    M.vecmatareisometries = @() false;
    
    
end

% Linear combination of tangent vectors
function d = lincomb(x, a1, d1, a2, d2) %#ok<INLSL>
    
    if nargin == 3
        d.L = a1*d1.L;
        d.R = a1*d1.R;
    elseif nargin == 5
        d.L = a1*d1.L + a2*d2.L;
        d.R = a1*d1.R + a2*d2.R;
    else
        error('Bad use of fixedrankfactory_2factors_subspace_projection.lincomb.');
    end
    
end

function A = uf(A)
    [L, unused, R] = svd(A, 0); %#ok
    A = L*R';
end

function omega = nested_sylvester(sym_mat, asym_mat)
    % omega=nested_sylvester(sym_mat,asym_mat)
    % This function solves the system of nested Sylvester equations:
    %
    %     X*sym_mat + sym_mat*X = asym_mat
    %     Omega*sym_mat+sym_mat*Omega = X
    % Mishra, Meyer, Bonnabel and Sepulchre, 'Fixed-rank matrix factorizations and Riemannian low-rank optimization'
    
    % Lses built-in lyap function, but does not exploit the fact that it's
    % twice the same sym_mat matrix that comes into play.
    
    X = lyap(sym_mat, -asym_mat);
    omega = lyap(sym_mat, -X);
    
end



