function M = obliquefactory(n, m, transposed)
% Returns a manifold struct to optimize over matrices w/ unit-norm columns.
%
% function M = obliquefactory(n, m)
% function M = obliquefactory(n, m, transposed)
%
% Oblique manifold: deals with matrices of size n x m such that each column
% has unit 2-norm, i.e., is a point on the unit sphere in R^n. The metric
% is such that the oblique manifold is a Riemannian submanifold of the
% space of nxm matrices with the usual trace inner product, i.e., the usual
% metric.
%
% If transposed is set to true (it is false by default), then the matrices
% are transposed: a point Y on the manifold is a matrix of size m x n and
% each row has unit 2-norm. It is the same geometry, just a different
% representation.
%
% See also: spherefactory

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 
%
%	July 16, 2013 (NB) :
%       Added 'transposed' option, mainly for ease of comparison with the
%       elliptope geometry.
%
%	Nov. 29, 2013 (NB) :
%       Added normalize_columns function to make it easier to exploit the
%       bsxfun formulation of column normalization, which avoids using for
%       loops and provides performance gains. The exponential still uses a
%       for loop.

    
    if ~exist('transposed', 'var') || isempty(transposed)
        transposed = false;
    end
    
    if transposed
        trnsp = @(X) X';
    else
        trnsp = @(X) X;
    end

    M.name = @() sprintf('Oblique manifold OB(%d, %d)', n, m);
    
    M.dim = @() (n-1)*m;
    
    M.inner = @(x, d1, d2) d1(:).'*d2(:);
    
    M.norm = @(x, d) norm(d(:));
    
    M.dist = @(x, y) norm(real(acos(sum(trnsp(x).*trnsp(y), 1))));
    
    M.typicaldist = @() pi*sqrt(m);
    
    M.proj = @(X, U) trnsp(projection(trnsp(X), trnsp(U)));
    
    M.tangent = M.proj;
    
    % For Riemannian submanifolds, converting a Euclidean gradient into a
    % Riemannian gradient amounts to an orthogonal projection.
	M.egrad2rgrad = M.proj;
    
    M.ehess2rhess = @ehess2rhess;
    function rhess = ehess2rhess(X, egrad, ehess, U)
        X = trnsp(X);
        egrad = trnsp(egrad);
        ehess = trnsp(ehess);
        U = trnsp(U);
        
        PXehess = projection(X, ehess);
        inners = sum(X.*egrad, 1);
        rhess = PXehess - bsxfun(@times, U, inners);
        
        rhess = trnsp(rhess);
    end
    
    M.exp = @exponential;
    % Exponential on the oblique manifold
    function y = exponential(x, d, t)
        x = trnsp(x);
        d = trnsp(d);
        
        if nargin < 3
            t = 1.0;
        end

        m = size(x, 2);
        y = zeros(size(x));
        if t ~= 0
            for i = 1 : m
                y(:, i) = sphere_exponential(x(:, i), d(:, i), t);
            end
        else
            y = x;
        end

        y = trnsp(y);
    end

    M.log = @logarithm;
    function v = logarithm(x1, x2)
        x1 = trnsp(x1);
        x2 = trnsp(x2);
        
        v = M.proj(x1, x2 - x1);
        dists = acos(sum(x1.*x2, 1));
        norms = sqrt(sum(v.^2, 1));
		factors = dists./norms;
		% factors(dists <= 1e-6) = 1;
		v = bsxfun(@times, v, factors);
        
        v = trnsp(v);
    end

    M.retr = @retraction;
    % Retraction on the oblique manifold
    function y = retraction(x, d, t)
        x = trnsp(x);
        d = trnsp(d);
        
        if nargin < 3
            t = 1.0;
        end

        m = size(x, 2);
        if t ~= 0
			y = normalize_columns(x + t*d);
        else
            y = x;
        end

        y = trnsp(y);
    end

    M.hash = @(x) ['z' hashmd5(x(:))];
    
    M.rand = @() trnsp(random(n, m));
    
    M.randvec = @(x) trnsp(randomvec(n, m, trnsp(x)));
    
    M.lincomb = @lincomb;
    
    M.zerovec = @(x) trnsp(zeros(n, m));
    
    M.transp = @(x1, x2, d) M.proj(x2, d);
    
    M.pairmean = @pairmean;
    function y = pairmean(x1, x2)
        y = trnsp(x1+x2);
        y = normalize_columns(y);
        y = trnsp(y);
    end

    % vec returns a vector representation of an input tangent vector which
    % is represented as a matrix. mat returns the original matrix
    % representation of the input vector representation of a tangent
    % vector. vec and mat are thus inverse of each other. They are
    % furthermore isometries between a subspace of R^nm and the tangent
    % space at x.
    vect = @(X) X(:);
    M.vec = @(x, u_mat) vect(trnsp(u_mat));
    M.mat = @(x, u_vec) trnsp(reshape(u_vec, [n, m]));
    M.vecmatareisometries = @() true;

end

% Given a matrix X, returns the same matrix but with each column scaled so
% that they have unit 2-norm.
function X = normalize_columns(X)
	norms = sqrt(sum(X.^2, 1));
	X = bsxfun(@times, X, 1./norms);
end

% Orthogonal projection of the ambient vector H onto the tangent space at X
function PXH = projection(X, H)

    % Compute the inner product between each vector H(:, i) with its root
    % point X(:, i), that is, X(:, i).' * H(:, i). Returns a row vector.
    inners = sum(X.*H, 1);
    
    % Subtract from H the components of the H(:, i)'s that are parallel to
    % the root points X(:, i).
    PXH = H - bsxfun(@times, X, inners);

    % % Equivalent but slow code:
    % m = size(X, 2);
    % PXH = zeros(size(H));
    % for i = 1 : m
    %     PXH(:, i) = H(:, i) - X(:, i) * (X(:, i)'*H(:, i));
    % end

end

% Exponential on the sphere.
function y = sphere_exponential(x, d, t)

    if nargin == 2
        t = 1.0;
    end
    
    td = t*d;
    
    nrm_td = norm(td);
    
    if nrm_td > 1e-6
        y = x*cos(nrm_td) + (td/nrm_td)*sin(nrm_td);
    else
        % if the step is too small, to avoid dividing by nrm_td, we choose
        % to approximate with this retraction-like step.
        y = x + td;
        y = y / norm(y);
    end

end

% Uniform random sampling on the sphere.
function x = random(n, m)

    x = normalize_columns(randn(n, m));

end

% Random normalized tangent vector at x.
function d = randomvec(n, m, x)

    d = randn(n, m);
    d = projection(x, d);
    d = d / norm(d(:));

end

% Linear combination of tangent vectors
function d = lincomb(x, a1, d1, a2, d2) %#ok<INUSL>

    if nargin == 3
        d = a1*d1;
    elseif nargin == 5
        d = a1*d1 + a2*d2;
    else
        error('Bad use of oblique.lincomb.');
    end

end
