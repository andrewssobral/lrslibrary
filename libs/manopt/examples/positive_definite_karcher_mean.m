function X = positive_definite_karcher_mean(A)
% Computes a Karcher mean of a collection of positive definite matrices.
%
% function X = positive_definite_karcher_mean(A)
%
% Input:  A 3D matrix A of size nxnxm such that each slice A(:,:,k) is a
%         positive definite matrix of size nxn.
% 
% Output: A positive definite matrix X of size nxn which is a Karcher mean
%         of the m matrices in A, that is, X minimizes the sum of squared
%         Riemannian distances to the matrices in A:
%            f(X) = sum_k=1^m .5*dist^2(X, A(:, :, k))
%         The distance is defined by the natural metric on the set of
%         positive definite matrices: dist(X,Y) = norm(logm(X\Y), 'fro').
% 
% This simple example is not the best way to compute Karcher means. Its
% purpose it to serve as base code to explore other algorithms. In
% particular, in the presence of large noise, this algorithm seems to not
% be able to reach points with a very small gradient norm. This may be
% caused by insufficient accuracy in the gradient computation.

% This file is part of Manopt and is copyrighted. See the license file.
% 
% Main author: Nicolas Boumal, Sept. 3, 2013
% Contributors:
% 
% Change log:
% 
    
    % Generate some random data to test the function if none is given.
    if ~exist('A', 'var') || isempty(A)
        n = 5;
        m = 10;
        A = zeros(n, n, m);
        ref = diag(max(.1, 1+.1*randn(n, 1)));
        for i = 1 : m
            noise = 0.01*randn(n);
            noise = (noise + noise')/2;
            [V D] = eig(ref + noise);
            A(:, :, i) = V*diag(max(.01, diag(D)))*V';
        end
    end
    
    % Retrieve the size of the problem:
    % There are m matrices of size nxn to average.
    n = size(A, 1);
    m = size(A, 3);
    assert(n == size(A, 2), ...
           ['The slices of A must be square, i.e., the ' ...
	        'first and second dimensions of A must be equal.']);
    
    % Our search space is the set of positive definite matrices of size n.
    % Notice that this is the only place we specify on which manifold we
    % wish to compute Karcher means. Replacing this factory for another
    % geometry will yield code to compute Karcher means on that other
    % manifold, provided that manifold is equipped with a dist function and
    % a logarithmic map log.
    M = sympositivedefinitefactory(n);
    
    % Define a problem structure, specifying the manifold M, the cost
    % function and its gradient.
    problem.M = M;
    problem.cost = @cost;
    problem.grad = @grad;
    
    % The functions below make many redundant computations. This
    % performance hit can be alleviated by using the caching system. We go
    % for a simple implementation here, as a tutorial example.
    
    % Cost function
    function f = cost(X)
        f = 0;
        for k = 1 : m
            f = f + M.dist(X, A(:, :, k))^2;
        end
        f = f/(2*m);
    end

    % Riemannian gradient of the cost function
    function g = grad(X)
        g = M.zerovec(X);
        for k = 1 : m
            % Update g in a linear combination of the form
            % g = g - [something]/m.
            g = M.lincomb(X, 1, g, -1/m, M.log(X, A(:, :, k)));
        end
    end
    
    % Execute some checks on the derivatives for early debugging.
    % These things can be commented out of course.
    % The slopes should agree on part of the plot at least. In this case,
    % it is sometimes necessary to inspect the plot visually to make the
    % call, but it is indeed correct.
    % checkgradient(problem);
    % pause;
    
    % Execute this if you want to force using a proper parallel vector
    % transport. This is not necessary. If you omit this, the default
    % vector transport is the identity map, which is (of course) cheaper
    % and seems to perform well in practice.
    % M.transp = M.paralleltransp;
    
    % Issue a call to a solver. Default options are selected.
    % Our initial guess is the first data point.
    X = trustregions(problem, A(:, :, 1));

end
