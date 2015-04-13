function [Z, P, X, A] = sparse_pca(A, m, gamma)
% Sparse principal component analysis based on optimization over Stiefel.
%
% [Z, P, X] = sparse_pca(A, m, gamma)
%
% We consider sparse PCA applied to a data matrix A of size pxn, where p is
% the number of samples (observations) and n is the number of variables
% (features). We attempt to extract m different components. The parameter
% gamma, which must lie between 0 and the largest 2-norm of a column of
% A, tunes the balance between best explanation of the variance of the data
% (gamma = 0, mostly corresponds to standard PCA) and best sparsity of the
% principal components Z (gamma maximal, Z is zero). The variables
% contained in the columns of A are assumed centered (zero-mean).
%
% The output Z of size nxm represents the principal components. There are m
% columns, each one of unit norm and capturing a prefered direction of the
% data, while trying to be sparse. P has the same size as Z and represents
% the sparsity pattern of Z. X is an orthonormal matrix of size pxm
% produced internally by the algorithm.
%
% With classical PCA, the variability captured by m components is
% sum(svds(A, m))
% With the outputted Z, which should be sparser than normal PCA, it is
% sum(svd(A*Z))
%
% The method is based on the maximization of a differentiable function over
% the Stiefel manifold of dimension pxm. Notice that this dimension is
% independent of n, making this method particularly suitable for problems
% with many variables but few samples (n much larger than p). The
% complexity of each iteration of the algorithm is linear in n as a result.
%
% The theory behind this code is available in the paper
% http://jmlr.org/papers/volume11/journee10a/journee10a.pdf
% Generalized Power Method for Sparse Principal Component Analysis, by
% Journee, Nesterov, Richtarik and Sepulchre, JMLR, 2010.
% This implementation is not equivalent to the one described in that paper
% (and is independent from their authors) but is close in spirit
% nonetheless. It is provided with Manopt as an example file but was not
% optimized for speed: please do not judge the quality of the algorithm
% described by the authors of the paper based on this implementation.

% This file is part of Manopt and is copyrighted. See the license file.
% 
% Main author: Nicolas Boumal, Dec. 24, 2013
% Contributors:
% 
% Change log:
% 

    % If no input is provided, generate random data for a quick demo
    if nargin == 0
        n = 100;
        p = 10;
        m = 2;

        % Data matrix
        A = randn(p, n);

        % Regularization parameter. This should be between 0 and the largest
        % 2-norm of a column of A.
        gamma = 1;
        
    elseif nargin ~= 3
        error('Please provide 3 inputs (or none for a demo).');
    end
    
    % Execute the main algorithm: it will compute a sparsity pattern P.
    [P, X] = sparse_pca_stiefel_l1(A, m, gamma);
    
    % Compute the principal components in accordance with the sparsity.
    Z = postprocess(A, P, X);

end


% Sparse PCA based on the block sparse PCA algorithm with l1-penalty as
% featured in the reference paper by Journee et al. This is not the same
% algorithm but it is the same cost function optimized over the same search
% space. We force N = eye(m).
function [P, X] = sparse_pca_stiefel_l1(A, m, gamma)
    
    [p, n] = size(A); %#ok<NASGU>

    % The optimization takes place over a Stiefel manifold whose dimension
    % is independent of n. This is especially useful when there are many
    % more variables than samples.
    St = stiefelfactory(p, m);
    problem.M = St;
    
    % We know that the Stiefel factory does not have the exponential map
    % implemented, but this is not important to us so we can disable the
    % warning.
    warning('off', 'manopt:stiefel:exp');

    % In this helper function, given a point 'X' on the manifold we check
    % whether the caching structure 'store' has been populated with
    % quantities that are useful to compute at X or not. If they were not,
    % then we compute and store them now.
    function store = prepare(X, store)
        if ~isfield(store, 'ready') || ~store.ready
            store.AtX = A'*X;
            store.absAtX = abs(store.AtX);
            store.pos = max(0, store.absAtX - gamma);
            store.ready = true;
        end
    end

    % Define the cost function here and set it in the problem structure.
    problem.cost = @cost;
    function [f store] = cost(X, store)
        store = prepare(X, store);
        pos = store.pos;
        f = -.5*norm(pos, 'fro')^2;
    end

    % Here, we chose to define the Euclidean gradient (egrad instead of
    % grad) : Manopt will take care of converting it to the Riemannian
    % gradient.
    problem.egrad = @egrad;
    function [G store] = egrad(X, store)
        if ~isfield(store, 'G')
            store = prepare(X, store);
            pos = store.pos;
            AtX = store.AtX;
            sgAtX = sign(AtX);
            factor = pos.*sgAtX;
            store.G = -A*factor;
        end
        G = store.G;
    end

    % checkgradient(problem);
    % pause;

    % The optimization happens here. To improve the method, it may be
    % interesting to investigate better-than-random initial iterates and,
    % possibly, to fine tune the parameters of the solver.
    X = trustregions(problem);

    % Compute the sparsity pattern by thresholding
    P = abs(A'*X) > gamma;
    
end


% This post-processing algorithm produces a matrix Z of size nxm matching
% the sparsity pattern P and representing sparse principal components for
% A. This is to be called with the output of the main algorithm. This
% algorithm is described in the reference paper by Journee et al.
function Z = postprocess(A, P, X)
    fprintf('Post-processing... ');
    counter = 0;
    maxiter = 1000;
    tolerance = 1e-8;
    while counter < maxiter
        Z = A'*X;
        Z(~P) = 0;
        Z = Z*diag(1./sqrt(diag(Z'*Z)));
        X = ufactor(A*Z);
        counter = counter + 1;
        if counter > 1 && norm(Z0-Z, 'fro') < tolerance*norm(Z0, 'fro')
            break;
        end
        Z0 = Z;
    end
    fprintf('done, in %d iterations (max = %d).\n', counter, maxiter);
end

% Returns the U-factor of the polar decomposition of X
function U = ufactor(X)
    [W S V] = svd(X, 0); %#ok<ASGLU>
    U = W*V';
end
