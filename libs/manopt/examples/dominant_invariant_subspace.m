function [X, info] = dominant_invariant_subspace(A, p)
% Returns an orthonormal basis of the dominant invariant p-subspace of A.
%
% function X = dominant_invariant_subspace(A, p)
%
% Input: A real, symmetric matrix A of size nxn and an integer p < n.
% Output: A real, orthonormal matrix X of size nxp such that trace(X'*A*X)
%         is maximized. That is, the columns of X form an orthonormal basis
%         of a dominant subspace of dimension p of A. These are thus
%         eigenvectors associated with the largest eigenvalues of A (in no
%         particular order). Sign is important: 2 is deemed a larger
%         eigenvalue than -5.
%
% The optimization is performed on the Grassmann manifold, since only the
% space spanned by the columns of X matters. The implementation is short to
% show how Manopt can be used to quickly obtain a prototype. To make the
% implementation more efficient, one might first try to use the caching
% system, that is, use the optional 'store' arguments in the cost, grad and
% hess functions. Furthermore, using egrad2rgrad and ehess2rhess is quick
% and easy, but not always efficient. Having a look at the formulas
% implemented in these functions can help rewrite the code without them,
% possibly more efficiently.

% This file is part of Manopt and is copyrighted. See the license file.
%
% Main author: Nicolas Boumal, July 5, 2013
% Contributors:
%
% Change log:
%
%   NB Dec. 6, 2013:
%       We specify a max and initial trust region radius in the options.
    
    % Generate some random data to test the function
    if ~exist('A', 'var') || isempty(A)
        A = randn(128);
        A = (A+A')/2;
    end
    if ~exist('p', 'var') || isempty(p)
        p = 3;
    end
    
    % Make sure the input matrix is square and symmetric
    n = size(A, 1);
	assert(isreal(A), 'A must be real.')
    assert(size(A, 2) == n, 'A must be square.');
    assert(norm(A-A', 'fro') < n*eps, 'A must be symmetric.');
	assert(p<=n, 'p must be smaller than n.');
    
    % Define the cost and its derivatives on the Grassmann manifold
    Gr = grassmannfactory(n, p);
    problem.M = Gr;
    problem.cost = @(X)    -trace(X'*A*X);
    problem.grad = @(X)    -2*Gr.egrad2rgrad(X, A*X);
    problem.hess = @(X, H) -2*Gr.ehess2rhess(X, A*X, A*H, H);
    
    % Execute some checks on the derivatives for early debugging.
    % These things can be commented out of course.
    % checkgradient(problem);
    % pause;
    % checkhessian(problem);
    % pause;
    
    % Issue a call to a solver. A random initial guess will be chosen and
    % default options are selected except for the ones we specify here.
    options.Delta_bar = 8*sqrt(p);
    [X, costX, info, options] = trustregions(problem, [], options); %#ok<ASGLU>
    
    fprintf('Options used:\n');
    disp(options);
    
    % For our information, Manopt can also compute the spectrum of the
    % Riemannian Hessian on the tangent space at (any) X. Computing the
    % spectrum at the solution gives us some idea of the conditioning of
    % the problem. If we were to implement a preconditioner for the
    % Hessian, this would also inform us on its performance.
    %
    % Notice that (typically) all eigenvalues of the Hessian at the
    % solution are positive, i.e., we find an isolated minimizer. If we
    % replace the Grassmann manifold by the Stiefel manifold, hence still
    % optimizing over orthonormal matrices but ignoring the invariance
    % cost(XQ) = cost(X) for all Q orthogonal, then we see
    % dim O(p) = p(p-1)/2 zero eigenvalues in the Hessian spectrum, making
    % the optimizer not isolated anymore.
    if Gr.dim() < 512
        evs = hessianspectrum(problem, X);
        stairs(sort(evs));
        title(['Eigenvalues of the Hessian of the cost function ' ...
               'at the solution']);
        xlabel('Eigenvalue number (sorted)');
        ylabel('Value of the eigenvalue');
    end

end
