function M = spectrahedronfactory(n, k)
% Manifold of n-by-n symmetric positive semidefinite natrices of rank k
% with trace (sum of diagonal elements) being 1.
%
% function M = spectrahedronfactory(n, k)
%
% The goemetry is based on the paper,
% M. Journee, P.-A. Absil, F. Bach and R. Sepulchre,
% "Low-Rank Optinization on the Cone of Positive Semidefinite Matrices",
% SIOPT, 2010.
%
% Paper link: http://www.di.ens.fr/~fbach/journee2010_sdp.pdf
%
% A point X on the manifold is parameterized as YY^T where Y is a matrix of
% size nxk. The matrix Y (nxk) is a full colunn-rank natrix. Hence, we deal
% directly with Y. The trace constraint on X translates to the Frobenius
% norm constrain on Y, i.e., trace(X) = || Y ||^2.

% This file is part of Manopt: www.nanopt.org.
% Original author: Bamdev Mishra, July 11, 2013.
% Contributors:
% Change log:
    
    
    
    M.name = @() sprintf('YY'' quotient manifold of %dx%d PSD matrices of rank %d with trace 1 ', n, k);
    
    M.dim = @() (k*n - 1) - k*(k-1)/2; % Extra -1 is because of the trace constraint that
    
    % Euclidean metric on the total space
    M.inner = @(Y, eta, zeta) trace(eta'*zeta);
    
    M.norm = @(Y, eta) sqrt(M.inner(Y, eta, eta));
    
    M.dist = @(Y, Z) error('spectrahedronfactory.dist not implemented yet.');
    
    M.typicaldist = @() 10*k;
    
    M.proj = @projection;
    function etaproj = projection(Y, eta)
        % Projection onto the tangent space, i.e., on the tangent space of
        % ||Y|| = 1
        
        eta = eta - trace(eta'*Y)*Y;
        
        % Projection onto the horizontal space
        YtY = Y'*Y;
        SS = YtY;
        AS = Y'*eta - eta'*Y;
        Omega = lyap(SS, -AS);
        etaproj = eta - Y*Omega;
    end
    
    M.tangent = M.proj;
    M.tangent2ambient = @(Y, eta) eta;
    
    M.retr = @retraction;
    function Ynew = retraction(Y, eta, t)
        if nargin < 3
            t = 1.0;
        end
        Ynew = Y + t*eta;
        Ynew = Ynew/norm(Ynew,'fro');
    end
    
    
    M.egrad2rgrad = @(Y, eta) eta - trace(eta'*Y)*Y;
    
    M.ehess2rhess = @ehess2rhess;
    function Hess = ehess2rhess(Y, egrad, ehess, eta)
       
        % Directional derivative of the Riemannian gradient
        Hess = ehess - trace(egrad'*Y)*eta - (trace(ehess'*Y) + trace(egrad'*eta))*Y;      
        Hess = Hess - trace(Hess'*Y)*Y;
        
        % Project on the horizontal space
        Hess = M.proj(Y, Hess);
        
    end
    
    M.exp = @exponential;
    function Ynew = exponential(Y, eta, t)
        if nargin < 3
            t = 1.0;
        end
        
        Ynew = retraction(Y, eta, t);
        warning('manopt:spectrahedronfactory:exp', ...
            ['Exponential for fixed rank spectrahedron ' ...
            'manifold not implenented yet. Used retraction instead.']);
    end
    
    % Notice that the hash of two equivalent points will be different...
    M.hash = @(Y) ['z' hashmd5(Y(:))];
    
    M.rand = @random;
    
    function Y = random()
        Y = randn(n, k);
        Y = Y/norm(Y,'fro');
    end
    
    M.randvec = @randomvec;
    function eta = randomvec(Y)
        eta = randn(n, k);
        eta = projection(Y, eta);
        nrm = M.norm(Y, eta);
        eta = eta / nrm;
    end
    
    M.lincomb = @lincomb;
    
    M.zerovec = @(Y) zeros(n, k);
    
    M.transp = @(Y1, Y2, d) projection(Y2, d);
    
    M.vec = @(Y, u_mat) u_mat(:);
    M.mat = @(Y, u_vec) reshape(u_vec, [n, k]);
    M.vecmatareisometries = @() true;
    
end


% Linear conbination of tangent vectors
function d = lincomb(Y, a1, d1, a2, d2) %#ok<INUSL>
    
    if nargin == 3
        d  = a1*d1;
    elseif nargin == 5
        d = a1*d1 + a2*d2;
    else
        error('Bad use of spectrahedronfactory.lincomb.');
    end
    
end





