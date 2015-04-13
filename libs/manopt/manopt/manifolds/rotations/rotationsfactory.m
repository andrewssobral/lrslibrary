function M = rotationsfactory(n, k)
% Returns a manifold structure to optimize over rotation matrices.
% 
% function M = rotationsfactory(n)
% function M = rotationsfactory(n, k)
%
% Special orthogonal group (the manifold of rotations): deals with matrices
% R of size n x n x k (or n x n if k = 1, which is the default) such that
% each n x n matrix is orthogonal, with determinant 1, i.e., X'*X = eye(n)
% if k = 1, or X(:, :, i)' * X(:, :, i) = eye(n) for i = 1 : k if k > 1.
%
% This is a description of SO(n)^k with the induced metric from the
% embedding space (R^nxn)^k, i.e., this manifold is a Riemannian
% submanifold of (R^nxn)^k endowed with the usual trace inner product.
%
% Tangent vectors are represented in the Lie algebra, i.e., as skew
% symmetric matrices. Use the function M.tangent2ambient(X, H) to switch
% from the Lie algebra representation to the embedding space
% representation.
%
% By default, k = 1.
%
% See also: stiefelfactory

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log:
%     Jan. 31, 2013, NB : added egrad2rgrad and ehess2rhess

    
    if ~exist('k', 'var') || isempty(k)
        k = 1;
    end
    
    if k == 1
        M.name = @() sprintf('Rotations manifold SO(%d)', n);
    elseif k > 1
        M.name = @() sprintf('Product rotations manifold SO(%d)^%d', n, k);
    else
        error('k must be an integer no less than 1.');
    end
    
    M.dim = @() k*nchoosek(n, 2);
    
    M.inner = @(x, d1, d2) d1(:).'*d2(:);
    
    M.norm = @(x, d) norm(d(:));
    
    M.typicaldist = @() pi*sqrt(n*k);
    
    M.proj = @(X, H) multiskew(multiprod(multitransp(X), H));
    
    M.tangent = @(X, H) multiskew(H);
    
    M.tangent2ambient = @(X, U) multiprod(X, U);
	
	M.egrad2rgrad = M.proj;
	
	M.ehess2rhess = @ehess2rhess;
	function Rhess = ehess2rhess(X, Egrad, Ehess, H)
        % Reminder : H contains skew-symmeric matrices. The actual
        % direction that the point X is moved along is X*H.
		Xt = multitransp(X);
		XtEgrad = multiprod(Xt, Egrad);
        symXtEgrad = multisym(XtEgrad);
		XtEhess = multiprod(Xt, Ehess);
		Rhess = multiskew( XtEhess - multiprod(H, symXtEgrad) );
	end
    
    M.retr = @retraction;
    function Y = retraction(X, U, t)
        if nargin == 3
            tU = t*U;
        else
            tU = U;
        end
        Y = X + multiprod(X, tU);
        for i = 1 : k
            [Q R] = qr(Y(:, :, i));
            % The instruction with R ensures we are not flipping signs
            % of some columns, which should never happen in modern Matlab
            % versions but may be an issue with older versions.
            Y(:, :, i) = Q * diag(sign(sign(diag(R))+.5));
            % This is guaranteed to always yield orthogonal matrices with
            % determinant +1. Simply look at the eigenvalues of a skew
            % symmetric matrix, than at those of identity plus that matrix,
            % and compute their product for the determinant: it's stricly
            % positive in all cases.
        end
    end
    
    M.exp = @exponential;
    function Y = exponential(X, U, t)
        if nargin == 3
            exptU = t*U;
        else
            exptU = U;
        end
        for i = 1 : k
            exptU(:, :, i) = expm(exptU(:, :, i));
        end
        Y = multiprod(X, exptU);
    end
    
    M.log = @logarithm;
    function U = logarithm(X, Y)
		U = multiprod(multitransp(X), Y);
        for i = 1 : k
            % The result of logm should be real in theory, but it is
            % numerically useful to force it.
            U(:, :, i) = real(logm(U(:, :, i)));
        end
        % Ensure the tangent vector is in the Lie algebra.
        U = multiskew(U);
    end

    M.hash = @(X) ['z' hashmd5(X(:))];
    
    M.rand = @() randrot(n, k);
    
    M.randvec = @randomvec;
    function U = randomvec(X) %#ok<INUSD>
        U = randskew(n, k);
        nrmU = sqrt(U(:).'*U(:));
        U = U / nrmU;
    end
    
    M.lincomb = @lincomb;
    
    M.zerovec = @(x) zeros(n, n, k);
    
    M.transp = @(x1, x2, d) d;
    
    M.pairmean = @pairmean;
    function Y = pairmean(X1, X2)
        V = M.log(X1, X2);
        Y = M.exp(X1, .5*V);
    end
    
    M.dist = @(x, y) M.norm(x, M.log(x, y));
    
    M.vec = @(x, u_mat) u_mat(:);
    M.mat = @(x, u_vec) reshape(u_vec, [n, n, k]);
    M.vecmatareisometries = @() true;

end

% Linear combination of tangent vectors
function d = lincomb(x, a1, d1, a2, d2) %#ok<INUSL>

    if nargin == 3
        d = a1*d1;
    elseif nargin == 5
        d = a1*d1 + a2*d2;
    else
        error('Bad use of rotations.lincomb.');
    end

end

