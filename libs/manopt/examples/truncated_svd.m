function [U S V info] = truncated_svd(A, p)
% Returns an SVD decomposition of A truncated to rank p.
%
% function [U S V] = truncated_svd(A, p)
%
% Input: A real matrix A of size mxn and an integer p <= min(m, n).
% Output: An orthonormal matrix U of size mxp, an orthonormal matrix Y of
%         size nxp and a diagonal matrix S of size pxp with nonnegative and
%         decreasing diagonal entries such that USV.' is the best rank p
%         approximation of A according to the Frobenius norm. All real.
%         This function produces an output akin to svds.
% 
% The decomposition is obtained by maximizing
%   f(U, V) = .5*norm(U'*A*V, 'fro')^2
% where U, V are orthonormal. Notice that f(U*Q, V*R) = f(U, V) for all
% Q, R orthogonal pxp matrices. Hence, only the column spaces of U and V
% matter and we may perform the optimization over a product of two
% Grassmannian manifolds.
%
% It is easy to show that maximizing f is equivalent to minimizing g with
%   g(U, V) = min_S norm(U*S*V' - A, 'fro')^2,
% which confirms that we are going for a best low-rank approximation of A.
% 
% The inner workings of the Grassmann manifold use the built-in svd
% function of Matlab but only for matrices of size mxp and nxp to
% re-orthonormalize them.
% 
% Notice that we are actually chasing a best fixed-rank approximation of a
% matrix, which is best obtained by working directly over a manifold of
% fixed-rank matrices. This is simply an example script to demonstrate some
% functionalities of the toolbox.
% 
% The code can be modified to accept a function handle for A(x) = A*x
% instead of a matrix A, which is often useful. This would further require
% a function handle At for the transpose of A, such that At(x) = A.'*x.

% This file is part of Manopt and is copyrighted. See the license file.
% 
% Main author: Nicolas Boumal, July 5, 2013
% Contributors:
% 
% Change log:
% 

    
    % Generate some random data to test the function if none is given.
    if ~exist('A', 'var') || isempty(A)
        A = randn(42, 60);
    end
    if ~exist('p', 'var') || isempty(p)
        p = 5;
    end
    
    % Retrieve the size of the problem and make sure the requested
    % approximation rank is at most the maximum possible rank.
    [m n] = size(A);
    assert(p <= min(m, n), 'p must be smaller than the smallest dimension of A.');
    
    % Define the cost and its derivatives on the Grassmann manifold
    tuple.U = grassmannfactory(m, p);
    tuple.V = grassmannfactory(n, p);
    % All of the code will work just as well if we ignore the invariance
    % property of the cost function indicated above and thus place U and V
    % on the Stiefel manifold (orthonormal matrices) instead of the
    % Grassmann manifold. Working on Stiefel is expected to be slower
    % though, partly because de search space is higher dimensional and
    % partly because the optimizers are not isolated.
    % tuple.U = stiefelfactory(m, p);
    % tuple.V = stiefelfactory(n, p);
    M = productmanifold(tuple);
    
    % Define a problem structure, specifying the manifold M, the cost
    % function and its derivatives. Here, to demonstrate the rapid
    % prototyping capabilities of Manopt, we directly define the Euclidean
    % gradient and the Euclidean Hessian egrad and ehess instead of the
    % Riemannian gradient and Hessian grad and hess. Manopt will take care
    % of the conversion. This automatic conversion is usually not
    % computationally optimal though, because much of the computations
    % involved in obtaining the gradient could be reused to obtain the
    % Hessian. After the prototyping stage, when efficiency becomes
    % important, it makes sense to define grad and hess rather than egrad
    % an ehess, and to use the caching system (the store structure).
    problem.M = M;
    problem.cost  = @cost;
    problem.egrad = @egrad;
    problem.ehess = @ehess;
    
    % The functions below make many redundant computations. This
    % performance hit can be alleviated by using the caching system.
    
    % Cost function
    function f = cost(X)
        U = X.U;
        V = X.V;
        f = -.5*norm(U'*A*V, 'fro')^2;
    end
    % Euclidean gradient of the cost function
    function g = egrad(X)
        U = X.U;
        V = X.V;
        AV = A*V;
        AtU = A'*U;
        g.U = -AV*(AV'*U);
        g.V = -AtU*(AtU'*V);
    end
    % Euclidean Hessian of the cost function
    function h = ehess(X, H)
        U = X.U;
        V = X.V;
        Udot = H.U;
        Vdot = H.V;
        AV = A*V;
        AtU = A'*U;
        AVdot = A*Vdot;
        AtUdot = A'*Udot;
        h.U = -(AVdot*AV'*U + AV*AVdot'*U + AV*AV'*Udot);
        h.V = -(AtUdot*AtU'*V + AtU*AtUdot'*V + AtU*AtU'*Vdot);
    end
    
    
    % Execute some checks on the derivatives for early debugging.
    % These things can be commented out of course.
    checkgradient(problem);
    pause;
    checkhessian(problem);
    pause;
    
    % Issue a call to a solver. A random initial guess will be chosen and
    % default options are selected. Here, we specify a maximum trust
    % region radius (which in turn induces an initial trust region radius).
    % Note that this is not required: default values are used if we omit
    % this. The diameter of the manifold scales like sqrt(2*p), hence the
    % form of our (empirical) choice.
    options.Delta_bar = 4*sqrt(2*p);
    [X Xcost info] = trustregions(problem, [], options); %#ok<ASGLU>
    U = X.U;
    V = X.V;
    
    % Finish the job by rotating U and V such that the middle matrix S can
    % be diagonal with nonnegative, decreasing entries. This requires a
    % small svd of size pxp.
    Spp = U'*A*V;
    [Upp Spp Vpp] = svd(Spp);
    U = U*Upp;
    S = Spp;
    V = V*Vpp;
    
    % For our information, Manopt can also compute the spectrum of the
    % Riemannian Hessian on the tangent space at (any) X. Computing the
    % spectrum at the solution gives us some idea of the conditioning of
    % the problem. If we were to implement a preconditioner for the
    % Hessian, this would also inform us on its performance.
    %
    % Notice that if the optimization is performed on a product of Stiefel
    % manifolds instead of a product of Grassmannians, the double
    % invariance under the orthogonal group O(p) will appear as twice
    % p*(p-1)/2, thus p*(p-1) zero eigenvalues in the spectrum of the
    % Hessian. This means that the minimizers are not isolated, which
    % typically hinders convergence of second order algorithms.
    if M.dim() < 512
        evs = hessianspectrum(problem, X);
        stairs(sort(evs));
        title(['Eigenvalues of the Hessian of the cost function ' ...
               'at the solution']);
    end

end
