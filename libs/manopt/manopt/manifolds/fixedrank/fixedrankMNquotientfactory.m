function M = fixedrankMNquotientfactory(m, n, k)
% Manifold of m-by-n matrices of rank k with quotient geometry.
%
% function M = fixedrankMNquotientfactory(m, n, k)
%
% This follows the quotient geometry described in the following paper:
% P.-A. Absil, L. Amodei and G. Meyer,
% "Two Newton methods on the manifold of fixed-rank matrices endowed
%  with Riemannian quotient geometries", arXiv, 2012.
%
% Paper link: http://arxiv.org/abs/1209.0068
%
% A point X on the manifold is represented as a structure with two
% fields: M and N. The matrix M (mxk) is orthonormal, while the matrix N
% (nxk) is full-rank.
%
% Tangent vectors are represented as a structure with two fields (M, N).

% This file is part of Manopt: www.manopt.org.
% Original author: Bamdev Mishra, Dec. 30, 2012.
% Contributors:
% Change log:


M.name = @() sprintf('MN'' quotient manifold of %dx%d matrices of rank %d', m, n, k);

M.dim = @() (m+n-k)*k;

% Choice of the metric is motivated by the symmetry present in the
% space
M.inner = @(X, eta, zeta) eta.M(:).'*zeta.M(:) + eta.N(:).'*zeta.N(:);

M.norm = @(X, eta) sqrt(M.inner(X, eta, eta));

M.dist = @(x, y) error('fixedrankMNquotientfactory.dist not implemented yet.');

M.typicaldist = @() 10*k;

symm = @(X) .5*(X+X');
stiefel_proj = @(M, H) H - M*symm(M'*H);

M.egrad2rgrad = @egrad2rgrad;
    function eta = egrad2rgrad(X, eta)
        eta.M = stiefel_proj(X.M, eta.M);
    end

M.ehess2rhess = @ehess2rhess;
    function Hess = ehess2rhess(X, egrad, ehess, eta)
        
        % Directional derivative of the Riemannian gradient
        Hess.M = ehess.M - eta.M*symm(X.M'*egrad.M);
        Hess.M = stiefel_proj(X.M, Hess.M);
        
        Hess.N = ehess.N;
        
        % Projection onto the horizontal space
        Hess = M.proj(X, Hess);
    end


M.proj = @projection;
    function etaproj = projection(X, eta)
        
        % Start by projecting the vector from Rmp x Rnp to the tangent
        % space to the total space, that is, eta.M should be in the
        % tangent space to Stiefel at X.M and eta.N is arbitrary.
        eta.M = stiefel_proj(X.M, eta.M);
        
        % Now project from the tangent space to the horizontal space, that
        % is, take care of the quotient.
        
        % First solve a Sylvester equation (A symm., B skew-symm.)
        A = X.N'*X.N + eye(k);
        B = eta.M'*X.M + eta.N'*X.N;
        B = B-B';
        omega = lyap(A, -B);
        
        % And project along the vertical space to the horizontal space.
        etaproj.M = eta.M + X.M*omega;
        etaproj.N = eta.N + X.N*omega;
        
    end

M.exp = @exponential;
    function Y = exponential(X, eta, t)
        if nargin < 3
            t = 1.0;
        end
        
        A = t*X.M'*eta.M;
        S = t^2*eta.M'*eta.M;
        Y.M = [X.M t*eta.M]*expm([A -S ; eye(k) A])*eye(2*k, k)*expm(-A);
        
        % re-orthonormalize (seems necessary from time to time)
        [Q R] = qr(Y.M, 0);
        Y.M = Q * diag(sign(diag(R)));
        
        Y.N = X.N + t*eta.N;
        
    end

% Factor M lives on the Stiefel manifold, hence we will reuse its
% random generator.
stiefelm = stiefelfactory(m, k);

M.retr = @retraction;
    function Y = retraction(X, eta, t)
        if nargin < 3
            t = 1.0;
        end
        
        Y.M = uf(X.M + t*eta.M); % This is a valid retraction
        Y.N = X.N + t*eta.N;   
    end

M.hash = @(X) ['z' hashmd5([X.M(:) ; X.N(:)])];

M.rand = @random;
    function X = random()
        X.M = stiefelm.rand();
        X.N = randn(n, k);
    end

M.randvec = @randomvec;
    function eta = randomvec(X)
        eta.M = randn(m, k);
        eta.N = randn(n, k);
        eta = projection(X, eta);
        nrm = M.norm(X, eta);
        eta.M = eta.M / nrm;
        eta.N = eta.N / nrm;
    end

M.lincomb = @lincomb;

M.zerovec = @(X) struct('M', zeros(m, k), 'N', zeros(n, k));

M.transp = @(x1, x2, d) projection(x2, d);

end


% Linear combination of tangent vectors
function d = lincomb(x, a1, d1, a2, d2) %#ok<INMSL>

if nargin == 3
    d.M = a1*d1.M;
    d.N = a1*d1.N;
elseif nargin == 5
    d.M = a1*d1.M + a2*d2.M;
    d.N = a1*d1.N + a2*d2.N;
else
    error('Bad use of fixedrankMNquotientfactory.lincomb.');
end

end


function A = uf(A)
[L, unused, R] = svd(A, 0);
A = L*R';
end