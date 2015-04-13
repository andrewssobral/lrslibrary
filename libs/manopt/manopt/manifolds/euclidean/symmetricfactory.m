function M = symmetricfactory(n, k)
% Returns a manifold struct to optimize over k symmetric matrices of size n
%
% function M = symmetricfactory(n)
% function M = symmetricfactory(n, k)
%
% Returns M, a structure describing the Euclidean space of n-by-n symmetric
% matrices equipped with the standard Frobenius distance and associated
% trace inner product, as a manifold for Manopt.
% By default, k = 1. If k > 1, points and vectors are stored in 3D matrices
% X of size nxnxk such that each slice X(:, :, i), for i = 1:k, is
% symmetric.

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Jan. 22, 2014.
% Contributors: 
% Change log: 
    
    if ~exist('k', 'var') || isempty(k)
        k = 1;
    end

    M.name = @() sprintf('(Symmetric matrices of size %d)^%d', n, k);
    
    M.dim = @() k*n*(n+1)/2;
    
    M.inner = @(x, d1, d2) d1(:).'*d2(:);
    
    M.norm = @(x, d) norm(d(:), 'fro');
    
    M.dist = @(x, y) norm(x(:)-y(:), 'fro');
    
    M.typicaldist = @() sqrt(k)*n;
    
    M.proj = @(x, d) multisym(d);
    
    M.egrad2rgrad = M.proj;
    
    %M.ehess2rhess = @(x, eg, eh, d) eh;
    
    M.tangent = @(x, d) d;
    
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
    
    M.rand = @() multisym(randn(n, n, k));
    
    M.randvec = @randvec;
    function u = randvec(x) %#ok<INUSD>
        u = multisym(randn(n, n, k));
        u = u / norm(u(:), 'fro');
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
    
    M.zerovec = @(x) zeros(n, n, k);
    
    M.transp = @(x1, x2, d) d;
    
    M.pairmean = @(x1, x2) .5*(x1+x2);
    
    M.vec = @(x, u_mat) u_mat(:);
    M.mat = @(x, u_vec) reshape(u_vec, [m, n]);
    M.vecmatareisometries = @() true;

end
