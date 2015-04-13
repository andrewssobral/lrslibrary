function M = complexcirclefactory(n)
% Returns a manifold struct to optimize over unit-modulus complex numbers.
%
% function M = complexcirclefactory()
% function M = complexcirclefactory(n)
%
% Description of vectors z in C^n (complex) such that each component z(i)
% has unit modulus. The manifold structure is the Riemannian submanifold
% structure from the embedding space R^2 x ... x R^2, i.e., the complex
% circle is identified with the unit circle in the real plane.
%
% By default, n = 1.
%
% See also spherecomplexfactory

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 
%
%   July 7, 2014 (NB): Added ehess2rhess function.
%
    
    if ~exist('n', 'var')
        n = 1;
    end

    M.name = @() sprintf('Complex circle (S^1)^%d', n);
    
    M.dim = @() n;
    
    M.inner = @(z, v, w) real(v'*w);
    
    M.norm = @(x, v) norm(v);
    
    M.dist = @(x, y) norm(acos(conj(x) .* y));
    
    M.typicaldist = @() pi*sqrt(n);
    
    M.proj = @(z, u) u - real( conj(u) .* z ) .* z;	
    
    M.tangent = M.proj;
    
    % For Riemannian submanifolds, converting a Euclidean gradient into a
    % Riemannian gradient amounts to an orthogonal projection.
	M.egrad2rgrad = M.proj;
    
    M.ehess2rhess = @ehess2rhess;
    function rhess = ehess2rhess(z, egrad, ehess, zdot)
        rhess = M.proj(z, ehess - real(z.*conj(egrad)).*zdot);
    end
    
    M.exp = @exponential;
    function y = exponential(z, v, t)
        if nargin <= 2
            t = 1.0;
        end

        y = zeros(n, 1);        
        tv = t*v;

        nrm_tv = abs(tv);
        
        % We need to distinguish between very small steps and the others.
        % For very small steps, we use a a limit version of the exponential
        % (which actually coincides with the retraction), so as to not
        % divide by very small numbers.
        mask = nrm_tv > 1e-6;
        y(mask) = z(mask).*cos(nrm_tv(mask)) + ...
                  tv(mask).*(sin(nrm_tv(mask))./nrm_tv(mask));
        y(~mask) = z(~mask) + tv(~mask);
        y(~mask) = y(~mask) ./ abs(y(~mask));
        
    end
    
    M.retr = @retraction;
    function y = retraction(z, v, t)
        if nargin <= 2
            t = 1.0;
        end
        y = z+t*v;
        y = y ./ abs(y);
    end

    M.log = @logarithm;
    function v = logarithm(x1, x2)
        v = M.proj(x1, x2 - x1);
        di = M.dist(x1, x2);
        nv = norm(v);
		v = v * (di / nv);
    end
    
    M.hash = @(z) ['z' hashmd5( [real(z(:)) ; imag(z(:))] ) ];
    
    M.rand = @random;
    function z = random()
        z = randn(n, 1) + 1i*randn(n, 1);
        z = z ./ abs(z);
    end
    
    M.randvec = @randomvec;
    function v = randomvec(z)
        % i*z(k) is a basis vector of the tangent vector to the k-th circle
        v = randn(n, 1) .* (1i*z);
        v = v / norm(v);
    end
    
    M.lincomb = @lincomb;
    
    M.zerovec = @(x) zeros(n, 1);
    
    M.transp = @(x1, x2, d) M.proj(x2, d);
    
    M.pairmean = @pairmean;
    function z = pairmean(z1, z2)
        z = z1+z2;
        z = z ./ abs(z);
    end

    M.vec = @(x, u_mat) [real(u_mat) ; imag(u_mat)];
    M.mat = @(x, u_vec) u_vec(1:n) + 1i*u_vec((n+1):end);
    M.vecmatareisometries = @() true;

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
