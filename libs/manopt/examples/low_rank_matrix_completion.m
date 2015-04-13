function low_rank_matrix_completion()
% Given partial observation of a low rank matrix, attempts to complete it.
%
% function low_rank_matrix_completion()
%
% This example demonstrates how to use the geometry factory for the
% embedded submanifold of fixed-rank matrices, fixedrankembeddedfactory.
% This geometry is described in the paper
% "Low-rank matrix completion by Riemannian optimization"
% Bart Vandereycken - SIAM Journal on Optimization, 2013.
%
% This can be a starting point for many optimization problems of the form:
%
% minimize f(X) such that rank(X) = k, size(X) = [m, n].
%
% Note that the code is long because it showcases quite a few features of
% Manopt: most of the code is optional.
%
% Input:  None. This example file generate random data.
% 
% Output: None.

% This file is part of Manopt and is copyrighted. See the license file.
% 
% Main author: Nicolas Boumal, July 15, 2014
% Contributors:
% 
% Change log:
% 
    
    % Random data generation. First, choose the size of the problem.
    % We will complete a matrix of size mxn of rank k:
    m = 200;
    n = 500;
    k = 10;
    % Generate a random mxn matrix A of rank k
    L = randn(m, k);
    R = randn(n, k);
    A = L*R';
    % Generate a random mask for observed entries: P(i, j) = 1 if the entry
    % (i, j) of A is observed, and 0 otherwise.
    fraction = 4 * k*(m+n-k)/(m*n);
    P = sparse(rand(m, n) <= fraction);
    % Hence, we know the nonzero entries in PA:
    PA = P.*A;
    
    
    
    
    
    
    
    % Pick the manifold of matrices of size mxn of fixed rank k.
    problem.M = fixedrankembeddedfactory(m, n, k);

    % Define the problem cost function. The input X is a structure with
    % fields U, S, V representing a rank k matrix as U*S*V'.
    % f(X) = 1/2 * || P.*(X-A) ||^2
    problem.cost = @cost;
    function f = cost(X)
        % Note that it is very much inefficient to explicitly construct the
        % matrix X in this way. Seen as we only need to know the entries
        % of Xmat corresponding to the mask P, it would be far more
        % efficient to compute those only.
        Xmat = X.U*X.S*X.V';
        f = .5*norm( P.*Xmat - PA , 'fro')^2;
    end

    % Define the Euclidean gradient of the cost function, that is, the
    % gradient of f(X) seen as a standard function of X.
    % nabla f(X) = P.*(X-A)
    problem.egrad = @egrad;
    function G = egrad(X)
        % Same comment here about Xmat.
        Xmat = X.U*X.S*X.V';
        G = P.*Xmat - PA;
    end

    % This is optional, but it's nice if you have it.
    % Define the Euclidean Hessian of the cost at X, along H, where H is
    % represented as a tangent vector: a structure with fields Up, Vp, M.
    % This is the directional derivative of nabla f(X) at X along Xdot:
    % nabla^2 f(X)[Xdot] = P.*Xdot
    problem.ehess = @euclidean_hessian;
    function ehess = euclidean_hessian(X, H)
        % The function tangent2ambient transforms H (a tangent vector) into
        % its equivalent ambient vector representation. The output is a
        % structure with fields U, S, V such that U*S*V' is an mxn matrix
        % corresponding to the tangent vector H. Note that there are no
        % additional guarantees about U, S and V. In particular, U and V
        % are not orthonormal.
        ambient_H = problem.M.tangent2ambient(X, H);
        Xdot = ambient_H.U*ambient_H.S*ambient_H.V';
        % Same comment here about explicitly constructing the ambient
        % vector as an mxn matrix Xdot: we only need its entries
        % corresponding to the mask P, and this could be computed
        % efficiently.
        ehess = P.*Xdot;
    end
    



    % Check consistency of the gradient and the Hessian. Useful if you
    % adapt this example for a new cost function and you would like to make
    % sure there is no mistake.
    % warning('off', 'manopt:fixedrankembeddedfactory:exp');
    % checkgradient(problem); pause;
    % checkhessian(problem); pause;
    
    
    
    
    
    % Compute an initial guess. Points on the manifold are represented as
    % structures with three fields: U, S and V. U and V need to be
    % orthonormal, S needs to be diagonal.
    [U, S, V] = svds(PA, k);
    X0.U = U;
    X0.S = S;
    X0.V = V;
    
    % Minimize the cost function using Riemannian trust-regions, starting
    % from the initial guess X0.
    X = trustregions(problem, X0);
    
    % The reconstructed matrix is X, represented as a structure with fields
    % U, S and V.
    Xmat = X.U*X.S*X.V';
    fprintf('||X-A||_F = %g\n', norm(Xmat - A, 'fro'));
    
    
    
    
    
    
    % Alternatively, we could decide to use a solver such as
    % steepestdescent or conjugategradient. These solvers need to solve a
    % line-search problem at each iteration. Standard line searches in
    % Manopt have generic purpose systems to do this. But for the problem
    % at hand, it so happens that we can rather accurately guess how far
    % the line-search should look, and it would be a waste to not use that.
    % Look up the paper referenced above for the mathematical explanation
    % of the code below.
    % 
    % To tell Manopt about this special information, we specify the
    % linesearch hint function in the problem structure. Notice that this
    % is not the same thing as specifying a linesearch function in the
    % options structure.
    % 
    % Both the SD and the CG solvers will detect that we
    % specify the hint function below, and they will use an appropriate
    % linesearch algorithm by default, as a result. Typically, they will
    % try the step t*H first, then if it does not satisfy an Armijo
    % criterion, they will decrease t geometrically until satisfaction or
    % failure.
    % 
    % Just like the cost, egrad and ehess functions, the linesearch
    % function could use a store structure if you like. The present code
    % does not use the store structure, which means quite a bit of the
    % computations are made redundantly, and as a result a better method
    % could appear slower. See the Manopt tutorial about caching when you
    % are ready to switch from a proof-of-concept code to an efficient
    % code.
    %
    % The inputs are X (a point on the manifold) and H, a tangent vector at
    % X that is assumed to be a descent direction. That is, there exists a
    % positive t such that f(Retraction_X(tH)) < f(X). The function below
    % is supposed to output a "t" that it is a good "guess" at such a t.
    problem.linesearch = @linesearch_helper;
    function t = linesearch_helper(X, H)
        % Note that you would not usually need the Hessian for this.
        residual_omega = nonzeros(problem.egrad(X));
        dir_omega      = nonzeros(problem.ehess(X, H));
        t = - dir_omega \ residual_omega ;
    end

    % Notice that for this solver, the Hessian is not needed.
    [Xcg, xcost, info, options] = conjugategradient(problem, X0);
    
    fprintf('Take a look at the options that CG used:\n');
    disp(options);
    fprintf('And see how many trials were made at each line search call:\n');
    info_ls = [info.linesearch];
    disp([info_ls.costevals]);

    
    
    
    
    
    
    
    fprintf('Try it again without the linesearch helper.\n');
    
    % Remove the linesearch helper from the problem structure.
    problem = rmfield(problem, 'linesearch');
    
    [Xcg, xcost, info, options] = conjugategradient(problem, X0);
    
    fprintf('Take a look at the options that CG used:\n');
    disp(options);
    fprintf('And see how many trials were made at each line search call:\n');
    info_ls = [info.linesearch];
    disp([info_ls.costevals]);
    
    
    
    
    
    

    % If the problem has a small enough dimension, we may (for analysis
    % purposes) compute the spectrum of the Hessian at a point X. This may
    % help in studying the conditioning of a problem. If you don't provide
    % the Hessian, Manopt will approximate the Hessian with finite
    % differences of the gradient and try to estimate its "spectrum" (it's
    % not a proper linear operator). This can give some intuition, but
    % should not be relied upon.
    if problem.M.dim() < 100
        fprintf('Computing the spectrum of the Hessian...');
        s = hessianspectrum(problem, X);
        hist(s);
    end
    
    
    
    
end
