function M = spherefactory(n, m)
% Returns a manifold struct to optimize over unit-norm vectors or matrices.
%
% function M = spherefactory(n)
% function M = spherefactory(n, m)
%
% Manifold of n-by-m real matrices of unit Frobenius norm.
% By default, m = 1, which corresponds to the unit sphere in R^n. The
% metric is such that the sphere is a Riemannian submanifold of the space
% of nxm matrices with the usual trace inner product, i.e., the usual
% metric.
% 
% See also: obliquefactory spherecomplexfactory

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 

    
    if ~exist('m', 'var')
        m = 1;
    end

    if m == 1
        M.name = @() sprintf('Sphere S^%d', n-1);
    else
        M.name = @() sprintf('Unit F-norm %dx%d matrices', n, m);
    end
    
    M.dim = @() n*m-1;
    
    M.inner = @(x, d1, d2) d1(:).'*d2(:);
    
    M.norm = @(x, d) norm(d, 'fro');
    
    M.dist = @(x, y) real(acos(x(:).'*y(:)));
    
    M.typicaldist = @() pi;
    
    M.proj = @(x, d) d - x*(x(:).'*d(:));
    
    M.tangent = M.proj;
	
    % For Riemannian submanifolds, converting a Euclidean gradient into a
    % Riemannian gradient amounts to an orthogonal projection.
	M.egrad2rgrad = M.proj;
	
	M.ehess2rhess = @ehess2rhess;
	function rhess = ehess2rhess(x, egrad, ehess, u)
        rhess = M.proj(x, ehess) - (x(:)'*egrad(:))*u;
	end
    
    M.exp = @exponential;
    
    M.retr = @retraction;

    M.log = @logarithm;
    function v = logarithm(x1, x2)
        v = M.proj(x1, x2 - x1);
        di = M.dist(x1, x2);
		nv = norm(v, 'fro');
		v = v * (di / nv);
    end
    
    M.hash = @(x) ['z' hashmd5(x(:))];
    
    M.rand = @() random(n, m);
    
    M.randvec = @(x) randomvec(n, m, x);
    
    M.lincomb = @lincomb;
    
    M.zerovec = @(x) zeros(n, m);
    
    M.transp = @(x1, x2, d) M.proj(x2, d);
    
    M.pairmean = @pairmean;
    function y = pairmean(x1, x2)
        y = x1+x2;
        y = y / norm(y, 'fro');
    end

    M.vec = @(x, u_mat) u_mat(:);
    M.mat = @(x, u_vec) reshape(u_vec, [n, m]);
    M.vecmatareisometries = @() true;

end

% Exponential on the sphere
function y = exponential(x, d, t)

    if nargin == 2
        t = 1;
    end
    
    td = t*d;
    
    nrm_td = norm(td, 'fro');
    
    if nrm_td > 1e-6
        y = x*cos(nrm_td) + td*(sin(nrm_td)/nrm_td);
    else
        % if the step is too small, to avoid dividing by nrm_td, we choose
        % to approximate with this retraction-like step.
        y = x + td;
        y = y / norm(y, 'fro');
    end

end

% Retraction on the sphere
function y = retraction(x, d, t)

    if nargin == 2
        t = 1;
    end
    
    y = x + t*d;
    y = y / norm(y, 'fro');

end

% Uniform random sampling on the sphere.
function x = random(n, m)

    x = randn(n, m);
    x = x/norm(x, 'fro');

end

% Random normalized tangent vector at x.
function d = randomvec(n, m, x)

    d = randn(n, m);
    d = d - x*(x(:).'*d(:));
    d = d / norm(d, 'fro');

end

% Linear combination of tangent vectors
function d = lincomb(x, a1, d1, a2, d2) %#ok<INUSL>

    if nargin == 3
        d = a1*d1;
    elseif nargin == 5
        d = a1*d1 + a2*d2;
    else
        error('Bad use of sphere.lincomb.');
    end

end
