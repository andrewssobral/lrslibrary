function M = fixedrankembeddedfactory(m, n, k)
% Manifold struct to optimize fixed-rank matrices w/ an embedded geometry.
%
% function M = fixedrankembeddedfactory(m, n, k)
%
% Manifold of m-by-n real matrices of fixed rank k. This follows the
% geometry described in this paper (which for now is the documentation):
% B. Vandereycken, "Low-rank matrix completion by Riemannian optimization",
% 2011.
%
% Paper link: http://arxiv.org/pdf/1209.3834.pdf
%
% A point X on the manifold is represented as a structure with three
% fields: U, S and V. The matrices U (mxk) and V (nxk) are orthonormal,
% while the matrix S (kxk) is any /diagonal/, full rank matrix.
% Following the SVD formalism, X = U*S*V'. Note that the diagonal entries
% of S are not constrained to be nonnegative.
%
% Tangent vectors are represented as a structure with three fields: Up, M
% and Vp. The matrices Up (mxk) and Vp (mxk) obey Up'*U = 0 and Vp'*V = 0.
% The matrix M (kxk) is arbitrary. Such a structure corresponds to the
% following tangent vector in the ambient space of mxn matrices:
%   Z = U*M*V' + Up*V' + U*Vp'
% where (U, S, V) is the current point and (Up, M, Vp) is the tangent
% vector at that point.
%
% Vectors in the ambient space are best represented as mxn matrices. If
% these are low-rank, they may also be represented as structures with
% U, S, V fields, such that Z = U*S*V'. Their are no resitrictions on what
% U, S and V are, as long as their product as indicated yields a real, mxn
% matrix.
%
% The chosen geometry yields a Riemannian submanifold of the embedding
% space R^(mxn) equipped with the usual trace (Frobenius) inner product.

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 
%
%	Feb. 20, 2014 (NB):
%       Added function tangent to work with checkgradient.
%   June 24, 2014 (NB):
%       A couple modifications following
%       Bart Vandereycken's feedback:
%       - The checksum (hash) was replaced for a faster alternative: it's a
%         bit less "safe" in that collisions could arise with higher
%         probability, but they're still very unlikely.
%       - The vector transport was changed.
%       The typical distance was also modified, hopefully giving the
%       trustregions method a better initial guess for the trust region
%       radius, but that should be tested for different cost functions too.
%    July 11, 2014 (NB):
%       Added ehess2rhess and tangent2ambient, supplied by Bart.
%    July 14, 2014 (NB):
%       Added vec, mat and vecmatareisometries so that hessianspectrum now
%       works with this geometry. Implemented the tangent function.
%       Made it clearer in the code and in the documentation in what format
%       ambient vectors may be supplied, and generalized some functions so
%       that they should now work with both accepted formats.
%       It is now clearly stated that for a point X represented as a
%       triplet (U, S, V), the matrix S needs to be diagonal.

    M.name = @() sprintf('Manifold of %dx%d matrices of rank %d', m, n, k);
    
    M.dim = @() (m+n-k)*k;
    
    M.inner = @(x, d1, d2) d1.M(:).'*d2.M(:) + d1.Up(:).'*d2.Up(:) ...
                                             + d1.Vp(:).'*d2.Vp(:);
    
    M.norm = @(x, d) sqrt(M.inner(x, d, d));
    
    M.dist = @(x, y) error('fixedrankembeddedfactory.dist not implemented yet.');
    
    M.typicaldist = @() M.dim();
    
    % Given Z in tangent vector format, projects the components Up and Vp
    % such that they satisfy the tangent space constraints up to numerical
    % errors. If Z was indeed a tangent vector at X, this should barely
    % affect Z (it would not at all if we had infinite numerical accuracy).
    M.tangent = @tangent;
    function Z = tangent(X, Z)
        Z.Up = Z.Up - X.U*(X.U'*Z.Up);
        Z.Vp = Z.Vp - X.V*(X.V'*Z.Vp);
    end

    % For a given ambient vector Z, applies it to a matrix W. If Z is given
    % as a matrix, this is straightfoward. If Z is given as a structure
    % with fields U, S, V such that Z = U*S*V', the product is executed
    % efficiently.
    function ZW = apply_ambient(Z, W)
        if ~isstruct(Z)
            ZW = Z*W;
        else
            ZW = Z.U*(Z.S*(Z.V'*W));
        end
    end

    % Same as apply_ambient, but applies Z' to W.
    function ZtW = apply_ambient_transpose(Z, W)
        if ~isstruct(Z)
            ZtW = Z'*W;
        else
            ZtW = Z.V*(Z.S'*(Z.U'*W));
        end
    end
    
    % Orthogonal projection of an ambient vector Z represented as an mxn
    % matrix or as a structure with fields U, S, V to the tangent space at
    % X, in a tangent vector structure format.
    M.proj = @projection;
    function Zproj = projection(X, Z)
            
        ZV = apply_ambient(Z, X.V);
        UtZV = X.U'*ZV;
        ZtU = apply_ambient_transpose(Z, X.U);

        Zproj.M = UtZV;
        Zproj.Up = ZV  - X.U*UtZV;
        Zproj.Vp = ZtU - X.V*UtZV';

    end

    M.egrad2rgrad = @projection;
    
    % Code supplied by Bart.
    % Given the Euclidean gradient at X and the Euclidean Hessian at X
    % along H, where egrad and ehess are vectors in the ambient space and H
    % is a tangent vector at X, returns the Riemannian Hessian at X along
    % H, which is a tangent vector.
    M.ehess2rhess = @ehess2rhess;
    function rhess = ehess2rhess(X, egrad, ehess, H)
        
        % Euclidean part
        rhess = projection(X, ehess);
        
        % Curvature part
        T = apply_ambient(egrad, H.Vp)/X.S;
        rhess.Up = rhess.Up + (T - X.U*(X.U'*T));
        T = apply_ambient_transpose(egrad, H.Up)/X.S;
        rhess.Vp = rhess.Vp + (T - X.V*(X.V'*T));
        
    end

    % Transforms a tangent vector Z represented as a structure (Up, M, Vp)
    % into a structure with fields (U, S, V) that represents that same
    % tangent vector in the ambient space of mxn matrices, as U*S*V'.
    % This matrix is equal to X.U*Z.M*X.V' + Z.Up*X.V' + X.U*Z.Vp'. The
    % latter is an mxn matrix, which could be too large to build
    % explicitly, and this is why we return a low-rank representation
    % instead. Note that there are no guarantees on U, S and V other than
    % that USV' is the desired matrix. In particular, U and V are not (in
    % general) orthonormal and S is not (in general) diagonal.
    % (In this implementation, S is identity, but this might change.)
    M.tangent2ambient = @tangent2ambient;
    function Zambient = tangent2ambient(X, Z)
        Zambient.U = [X.U*Z.M + Z.Up, X.U];
        Zambient.S = eye(2*k);
        Zambient.V = [X.V, Z.Vp];
    end
    
    % This retraction is second order, following general results from
    % Absil, Malick, "Projection-like retractions on matrix manifolds",
    % SIAM J. Optim., 22 (2012), pp. 135-158.
    M.retr = @retraction;
    function Y = retraction(X, Z, t)
        if nargin < 3
            t = 1.0;
        end

        % See personal notes June 28, 2012 (NB)
        [Qu, Ru] = qr(Z.Up, 0);
        [Qv, Rv] = qr(Z.Vp, 0);
        
        % Calling svds or svd should yield the same result, but BV
        % advocated svd is more robust, and it doesn't change the
        % asymptotic complexity to call svd then trim rather than call
        % svds. Also, apparently Matlab calls ARPACK in a suboptimal way
        % for svds in this scenario.
        % [Ut St Vt] = svds([X.S+t*Z.M , t*Rv' ; t*Ru , zeros(k)], k);
        [Ut, St, Vt] = svd([X.S+t*Z.M , t*Rv' ; t*Ru , zeros(k)]);
        
        Y.U = [X.U Qu]*Ut(:, 1:k);
        Y.V = [X.V Qv]*Vt(:, 1:k);
        Y.S = St(1:k, 1:k) + eps*eye(k);
        
        % equivalent but very slow code
        % [U S V] = svds(X.U*X.S*X.V' + t*(X.U*Z.M*X.V' + Z.Up*X.V' + X.U*Z.Vp'), k);
        % Y.U = U; Y.V = V; Y.S = S;
        
    end
    
    M.exp = @exponential;
    function Y = exponential(X, Z, t)
        if nargin < 3
            t = 1.0;
        end
        Y = retraction(X, Z, t);
        warning('manopt:fixedrankembeddedfactory:exp', ...
               ['Exponential for fixed rank ' ...
                'manifold not implemented yet. Used retraction instead.']);
    end

    % Less safe but much faster checksum, June 24, 2014.
    % Older version right below.
    M.hash = @(X) ['z' hashmd5([sum(X.U(:)) ; sum(X.S(:)); sum(X.V(:)) ])];
    %M.hash = @(X) ['z' hashmd5([X.U(:) ; X.S(:) ; X.V(:)])];
    
    M.rand = @random;
    % Factors U and V live on Stiefel manifolds, hence we will reuse
    % their random generator.
    stiefelm = stiefelfactory(m, k);
    stiefeln = stiefelfactory(n, k);
    function X = random()
        X.U = stiefelm.rand();
        X.V = stiefeln.rand();
        X.S = diag(sort(rand(k, 1), 1, 'descend'));
    end
    
    % Generate a random tangent vector at X.
    % TODO: consider a possible imbalance between the three components Up,
    % Vp and M, when m, n and k are widely different (which is typical).
    M.randvec = @randomvec;
    function Z = randomvec(X)
        Z.Up = randn(m, k);
        Z.Vp = randn(n, k);
        Z.M  = randn(k);
        Z = tangent(X, Z);
        nrm = M.norm(X, Z);
        Z.Up = Z.Up / nrm;
        Z.Vp = Z.Vp / nrm;
        Z.M  = Z.M  / nrm;
    end
    
    M.lincomb = @lincomb;
    
    M.zerovec = @(X) struct('Up', zeros(m, k), 'M', zeros(k, k), ...
                                                        'Vp', zeros(n, k));
    
    % New vector transport on June 24, 2014 (as indicated by Bart)
    % Reference: Absil, Mahony, Sepulchre 2008 section 8.1.3:
    % For Riemannian submanifolds of a Euclidean space, it is acceptable to
    % transport simply by orthogonal projection of the tangent vector
    % translated in the ambient space.
    M.transp = @project_tangent;
    function Z2 = project_tangent(X1, X2, Z1)
        Z2 = projection(X2, tangent2ambient(X1, Z1));
    end


    M.vec = @vec;
    function Zvec = vec(X, Z)
        Zamb = tangent2ambient(X, Z);
        Zamb_mat = Zamb.U*Zamb.S*Zamb.V';
        Zvec = Zamb_mat(:);
    end
    M.mat = @(X, Zvec) projection(X, reshape(Zvec, [m, n]));
    M.vecmatareisometries = @() true;

end

% Linear combination of tangent vectors
function d = lincomb(x, a1, d1, a2, d2) %#ok<INUSL>

    if nargin == 3
        d.Up = a1*d1.Up;
        d.Vp = a1*d1.Vp;
        d.M  = a1*d1.M;
    elseif nargin == 5
        d.Up = a1*d1.Up + a2*d2.Up;
        d.Vp = a1*d1.Vp + a2*d2.Vp;
        d.M  = a1*d1.M  + a2*d2.M;
    else
        error('fixedrank.lincomb takes either 3 or 5 inputs.');
    end

end
