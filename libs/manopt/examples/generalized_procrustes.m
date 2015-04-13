function [A R] = generalized_procrustes(A_measure)
% Rotationally align clouds of points (generalized Procrustes problem)
%
% function X = generalized_procrustes(A_measure)
%
% The input is a 3D matrix A_measure of size nxmxN. Each of the N slices
% A_measure(:, :, i) is a cloud of m points in R^n. These clouds are
% assumed to be (noisy) rotated versions of a reference cloud Atrue.
% This algorithm tries to find the optimal rotations to apply to the
% individual clouds such that they will match each other as much as
% possible following a least-squares cost.
%
% The output A is an estimate of the cloud Atrue (up to rotation). The
% output R is a 3D matrix of size nxnxN containing the rotation matrices
% such that R(:, :, i) * A is approximately equal to A_measure(:, :, i).

% This file is part of Manopt and is copyrighted. See the license file.
%
% Main author: Nicolas Boumal, July 8, 2013
% Contributors:
%
% Change log:
%   
    
    if ~exist('A_measure', 'var')
        % Generate random data to test the method.
        % There are N clouds of m points in R^n. Each of them is a noisy,
        % rotated version of a reference cloud A. Rotations are uniformly
        % random and noise on each rotated cloud is iid normal with
        % standard deviation sigma.
        n = 3;
        m = 10;
        N = 50;
        % The reference cloud
        Atrue = randn(n, m);
        % A 3D meatrix containing the N measured clouds
        sigma = .3;
        A_measure = multiprod(randrot(n, N), Atrue) + sigma*randn(n, m, N);
    else
        [n, m, N] = size(A_measure);
    end
    
    % Construct a manifold structure representing the product of groups of
    % rotations with the Euclidean space for A. We optimize simultaneously
    % for the reference cloud and for the rotations that affect each of the
    % measured clouds. Notice that there is a group invariance because
    % there is no way of telling which orientation the reference cloud
    % should be in.
    tuple.R = rotationsfactory(n, N);
    tuple.A = euclideanfactory(n, m);
    M = productmanifold(tuple);

    % Define the cost function here. Points on the manifold M are
    % structures with fields X.A and X.R, containing matrices of sizes
    % respectively nxm and nxnxN. The store structure (the caching system)
    % is used to keep the residue matrix E in memory, as it is also used in
    % the computation of the gradient and of the Hessian. This way, we
    % prevent redundant computations.
    function [f store] = cost(X, store)
        if ~isfield(store, 'E')
            R = X.R;
            A = X.A;
            store.E = multiprod(R, A) - A_measure;
        end
        E = store.E;
        f = (E(:)'*E(:))/(2*N);
    end

    % Riemannian gradient of the cost function.
    function [g store] = grad(X, store)
        R = X.R;
        A = X.A;
        if ~isfield(store, 'E')
            [~, store] = cost(X, store);
        end
        E = store.E;
        % Compute the Euclidean gradient of the cost wrt the rotations R
        % and wrt the cloud A,
        egrad.R = multiprod(E, A'/N);
        egrad.A = A - mean(multiprod(multitransp(R), A_measure), 3);
        % then transform this Euclidean gradient into the Riemannian
        % gradient.
        g = M.egrad2rgrad(X, egrad);
        store.egrad = egrad;
    end

    % It is not necessary to define the Hessian of the cost. We do it
    % mostly to illustrate how to do it and to study the spectrum of the
    % Hessian at the solution (see further down).
    function [h store] = hess(X, Xdot, store)
        R = X.R;
        A = X.A;
        % Careful: tangent vectors on the rotation group are represented as
        % skew symmetric matrices. To obtain the corresponding vectors in
        % the ambient space, we need a little transformation. This
        % transformation is typically not needed when we compute the
        % formulas for the gradient and the Hessian directly in Riemannian
        % form instead of resorting the egrad2rgrad and ehess2rhess. These
        % latter tools are convenient for prototyping but are not always
        % the most efficient form to execute the computations.
        Rdot = tuple.R.tangent2ambient(R, Xdot.R);
        Adot = Xdot.A;
        if ~isfield(store, 'egrad')
            [~, store] = grad(X, store);
        end
        E = store.E;
        egrad = store.egrad;
        
        ehess.R = multiprod(multiprod(Rdot, A) + multiprod(R, Adot), A') + ...
                  multiprod(E, Adot');
        ehess.R = ehess.R / N;
        ehess.A = Adot-mean(multiprod(multitransp(Rdot), A_measure), 3);
        
        h = M.ehess2rhess(X, egrad, ehess, Xdot);
    end

    % Setup the problem structure with manifold M and cost+grad functions.
    problem.M = M;
    problem.cost = @cost;
    problem.grad = @grad;
    problem.hess = @hess;

    % For debugging, it's always nice to check the gradient a few times.
    % checkgradient(problem);
    % pause;
    % checkhessian(problem);
    % pause;
    
    % Call a solver on our problem. This can probably be much improved if a
	% clever initial guess is used instead of a random one.
    X = trustregions(problem);
    A = X.A;
    R = X.R;
    
    % To evaluate the performance of the algorithm, see how well Atrue (the
    % reference cloud) matches A (the found cloud). Since the recovery is
    % up to rotation, apply Kabsch algorithm (or standard Procrustes),
    % i.e., compute the polar factorization to best align Atrue and A.
    if exist('Atrue', 'var')
        [U, ~, V] = svd(Atrue*A');
        Ahat = (U*V')*A;
        fprintf('Registration error: %g.\n', norm(Atrue-Ahat, 'fro'));
    end
    
    % Plot the spectrum of the Hessian at the solution found.
    % Notice that the invariance of f under a rotation yields dim SO(n),
    % that is, n*(n-1)/2 zero eigenvalues in the Hessian spectrum at the
    % solution. This indicates that critical points are not isolated and
    % can theoretically prevent quadratic convergence. One solution to
    % circumvent this would be to fix one rotation arbitrarily. Another
    % solution would be to work on a quotient manifold. Both can be
    % achieved in Manopt: they simply require a little more work on the
    % manifold description side.
    if M.dim() <= 512
        stairs(sort(hessianspectrum(problem, X)));
        title('Spectrum of the Hessian at the solution found.');
        xlabel('Eigenvalue number (sorted)');
        ylabel('Value of the eigenvalue');
    end
    
end
