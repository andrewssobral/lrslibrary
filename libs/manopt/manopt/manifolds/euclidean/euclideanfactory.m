function M = euclideanfactory(m, n)
% Returns a manifold struct to optimize over m-by-n matrices.
%
% function M = euclideanfactory(m, n)
%
% Returns M, a structure describing the Euclidean space of m-by-n matrices
% equipped with the standard Frobenius distance and associated trace inner
% product as a manifold for Manopt.

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 
%  July 5, 2013 (NB): added egred2rgrad, ehess2rhess, mat, vec, tangent.

    
    if ~exist('n', 'var') || isempty(n)
        n = 1;
    end

    M.name = @() sprintf('Euclidean space R^(%dx%d)', m, n);
    
    M.dim = @() m*n;
    
    M.inner = @(x, d1, d2) d1(:).'*d2(:);
    
    M.norm = @(x, d) norm(d, 'fro');
    
    M.dist = @(x, y) norm(x-y, 'fro');
    
    M.typicaldist = @() sqrt(m*n);
    
    M.proj = @(x, d) d;
    
    M.egrad2rgrad = @(x, g) g;
    
    M.ehess2rhess = @(x, eg, eh, d) eh;
    
    M.tangent = M.proj;
    
    M.exp = @exp;
    function y = exp(x, d, t)
        if nargin == 3
            y = x + t*d;
        else
            y = x + d;
        end
    end
    
    M.retr = M.exp;
	
	M.log = @(x, y) y-x;

    M.hash = @(x) ['z' hashmd5(x(:))];
    
    M.rand = @() randn(m, n);
    
    M.randvec = @randvec;
    function u = randvec(x) %#ok<INUSD>
        u = randn(m, n);
        u = u / norm(u, 'fro');
    end
    
    M.lincomb = @lincomb;
    function v = lincomb(x, a1, d1, a2, d2) %#ok<INUSL>
        if nargin == 3
            v = a1*d1;
        elseif nargin == 5
            v = a1*d1 + a2*d2;
        else
            error('Bad usage of euclidean.lincomb');
        end
    end
    
    M.zerovec = @(x) zeros(m, n);
    
    M.transp = @(x1, x2, d) d;
    
    M.pairmean = @(x1, x2) .5*(x1+x2);
    
    M.vec = @(x, u_mat) u_mat(:);
    M.mat = @(x, u_vec) reshape(u_vec, [m, n]);
    M.vecmatareisometries = @() true;

end
