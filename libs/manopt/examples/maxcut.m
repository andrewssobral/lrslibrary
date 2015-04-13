function [x cutvalue cutvalue_upperbound Y] = maxcut(L, r)
% Algorithm to (try to) compute a maximum cut of a graph, via SDP approach.
% 
% function x = maxcut(L)
% function [x cutvalue cutvalue_upperbound Y] = maxcut(L, r)
%
% L is the Laplacian matrix describing the graph to cut. The Laplacian of a
% graph is L = D - A, where D is the diagonal degree matrix (D(i, i) is the
% sum of the weights of the edges adjacent to node i) and A is the
% symmetric adjacency matrix of the graph (A(i, j) = A(j, i) is the weight
% of the edge joining nodes i and j). If L is sparse, this will be
% exploited.
%
% If the graph has n nodes, then L is nxn and the output x is a vector of
% length n such that x(i) is +1 or -1. This partitions the nodes of the
% graph in two classes, in an attempt to maximize the sum of the weights of
% the edges that go from one class to the other (MAX CUT problem).
%
% cutvalue is the sum of the weights of the edges 'cut' by the partition x.
%
% If the algorithm reached the global optimum of the underlying SDP
% problem, then it produces an upperbound on the maximum cut value. This
% value is returned in cutvalue_upperbound if it is found. Otherwise, that
% output is set to NaN.
%
% If r is specified (by default, r = n), the algorithm will stop at rank r.
% This may prevent the algorithm from reaching a globally optimal solution
% for the underlying SDP problem (but can greatly help in keeping the
% execution time under control). If a global optimum of the SDP is reached
% before rank r, the algorithm will stop of course.
%
% Y is a matrix of size nxp, with p <= r, such that X = Y*Y' is the best
% solution found for the underlying SDP problem. If cutvalue_upperbound is
% not NaN, then Y*Y' is optimal for the SDP and cutvalue_upperbound is its
% cut value.
% 
% By Goemans and Williamson 1995, it is known that if the optimal value of
% the SDP is reached, then the returned cut, in expectation, is at most at
% a fraction 0.878 of the optimal cut. (This is not exactly valid because
% we do not use random projection here; sign(Y*randn(size(Y, 2), 1)) will
% give a cut that respects this statement -- it's usually worse though).
%
% The algorithm is essentially that of:
% Journee, Bach, Absil and Sepulchre, 2010
% Low-rank optimization on the code of positive semidefinite matrices.
%
% It is itself based on the famous SDP relaxation of MAX CUT:
% Goemans and Williamson, 1995
% Improved approximation algorithms for maximum cut and satisfiability
% problems using semidefinite programming.

% This file is part of Manopt and is copyrighted. See the license file.
%
% Main author: Nicolas Boumal, July 18, 2013
% Contributors:
%
% Change log:
%   


    % If no inputs are provided, generate a random Laplacian.
    % This is for illustration purposes only.
    if ~exist('L', 'var') || isempty(L)
        n = 20;
        A = triu(randn(n) <= .4, 1);
        A = A+A';
        D = diag(sum(A, 2));
        L = D-A;
    end


    n = size(L, 1);
    assert(size(L, 2) == n, 'L must be square.');

    if ~exist('r', 'var') || isempty(r) || r > n
        r = n;
    end
    
    % We will let the rank increase. Each rank value will generate a cut.
    % We have to go up in the rank to eventually find a certificate of SDP
    % optimality. This in turn will give us an upperbound on the MAX CUT
    % value and assure us that we're doing well, according to Goemans and
    % Williamson's argument. In practice though, the good cuts often come
    % up for low rank values, so we better keep track of the best one.
    best_x = ones(n, 1);
    best_cutvalue = 0;
    cutvalue_upperbound = NaN;
    
    time = [];
    cost = [];
    
    for rr = 2 : r
        
        manifold = elliptopefactory(n, rr);
        
        if rr == 2
            
            % At first, for rank 2, generate a random point.
            Y0 = manifold.rand();
             
        else
            
            % To increase the rank, we could just add a column of zeros to
            % the Y matrix. Unfortunately, this lands us in a saddle point.
            % To escape from the saddle, we may compute an eigenvector of
            % Sy associated to a negative eigenvalue: that will yield a
            % (second order) descent direction Z. See Journee et al ; Sy is
            % linked to dual certificates for the SDP.
            Y0 = [Y zeros(n, 1)];
            LY0 = L*Y0;
            Dy = spdiags(sum(LY0.*Y0, 2), 0, n, n);
            Sy = (Dy - L)/4;
            % Find the smallest (the "most negative") eigenvalue of Sy.
            [v, s] = eigs(Sy, 1, 'SA');
            % If there is no negative eigenvalue for Sy, than we are not at
            % a saddle point: we're actually done!
            if s >= -1e-8
                % We can stop here: we found the global optimum of the SDP,
                % and hence the reached cost is a valid upper bound on the
                % maximum cut value.
                cutvalue_upperbound = max(-[info.cost]);
                break;
            end
            
            % This is our escape direction.
            Z = manifold.proj(Y0, [zeros(n, rr-1) v]);
            
            % % These instructions can be uncommented to see what the cost
            % % function looks like at a saddle point. But will require the
            % % problem structure which is not defined here: see the helper
            % % function.
            % plotprofile(problem, Y0, Z, linspace(-1, 1, 101));
            % drawnow; pause;
            
            % Now make a step in the Z direction to escape from the saddle.
            % It is not obvious that it is ok to do a unit step ... perhaps
            % need to be cautious here with the stepsize. It's not too
            % critical though: the important point is to leave the saddle
            % point. But it's nice to guarantee monotone decrease of the
            % cost, and we can't do that with a constant step (at least,
            % not without a proper argument to back it up).
            stepsize = 1;
            Y0 = manifold.retr(Y0, Z, stepsize);
            
        end
        
        % Use the Riemannian optimization based algorithm lower in this
        % file to reach a critical point (typically a local optimizer) of
        % the max cut cost with fixed rank, starting from Y0.
        [Y info] = maxcut_fixedrank(L, Y0);
        
        % Some info logging.
        thistime = [info.time];
        if ~isempty(time)
            thistime = time(end) + thistime;
        end
        time = [time thistime]; %#ok<AGROW>
        cost = [cost [info.cost]]; %#ok<AGROW>

        % Time to turn the matrix Y into a cut.
        % We can either do the random rounding as follows:
        % x = sign(Y*randn(rr, 1));
        % or extract the "PCA direction" of the points in Y and cut
        % orthogonally to that direction, as follows:
        [u, ~, ~] = svds(Y, 1);
        x = sign(u);

        cutvalue = (x'*L*x)/4;
        if cutvalue > best_cutvalue
            best_x = x;
            best_cutvalue = cutvalue;
        end
        
    end
    
    x = best_x;
    cutvalue = best_cutvalue;
    
    plot(time, -cost, '.-');
    xlabel('Time [s]');
    ylabel('Relaxed cut value');
    title('The relaxed cut value is an upper bound on the optimal cut value.');

end


function [Y info] = maxcut_fixedrank(L, Y)
% Try to solve the (fixed) rank r relaxed max cut program, based on the
% Laplacian of the graph L and an initial guess Y. L is nxn and Y is nxr.

    [n r] = size(Y);
    assert(all(size(L) == n));
    
    % The fixed rank elliptope geometry describes symmetric, positive
    % semidefinite matrices of size n with rank r and all diagonal entries
    % are 1.
    manifold = elliptopefactory(n, r);
    
    % % If you want to compare the performance of the elliptope geometry
    % % against the (conceptually simpler) oblique manifold geometry,
    % % uncomment this line.
    % manifold = obliquefactory(r, n, true);
    
    problem.M = manifold;
    
    % % For rapid prototyping, these lines suffice to describe the cost
    % % function and its gradient and Hessian (here expressed using the
    % % Euclidean gradient and Hessian).
    % problem.cost = @(Y)  -trace(Y'*L*Y)/4;
    % problem.egrad = @(Y) -(L*Y)/2;
    % problem.ehess = @(Y, U) -(L*U)/2;
    
    % Instead of the prototyping version, the functions below describe the
    % cost, gradient and Hessian using the caching system (the store
    % structure). This alows to execute exactly the required number of
    % multiplications with the matrix L. These multiplications are counted
    % using the Lproducts_counter and registered for each iteration in the
    % info structure outputted by solvers, via the statsfun function.
    % Notice that we do not use the store structure to count: this does not
    % behave well in general and is not advised.
    
    Lproducts_counter = 0;

    % For every visited point Y, we will need L*Y. This function makes sure
    % the quantity L*Y is available, but only computes it if it wasn't
    % already computed.
    function store = prepare(Y, store)
        if ~isfield(store, 'LY')
            store.LY = L*Y;
            Lproducts_counter = Lproducts_counter + 1;
        end
    end

    problem.cost = @cost;
    function [f store] = cost(Y, store)
        store = prepare(Y, store);
        LY = store.LY;
        f = -(Y(:)'*LY(:))/4; % = -trace(Y'*LY)/4;
    end

    problem.grad = @grad;
    function [g store] = grad(Y, store)
        store = prepare(Y, store);
        LY = store.LY;
        g = manifold.egrad2rgrad(Y, -LY/2);
    end

    problem.hess = @hess;
    function [h store] = hess(Y, U, store)
        store = prepare(Y, store);
        LY = store.LY;
        LU = L*U;
        Lproducts_counter = Lproducts_counter + 1;
        h = manifold.ehess2rhess(Y, -LY/2, -LU/2, U);
    end

    % statsfun is called exactly once after each iteration (including after
    % the evaluation of the cost at the initial guess). We then register
    % the value of the Lproducts counter (which counts how many product
    % were needed since the last iteration), and reset it to zero.
    options.statsfun = @statsfun;
    function stats = statsfun(problem, Y, stats, store) %#ok
        stats.Lproducts = Lproducts_counter;
        Lproducts_counter = 0;
    end
    

    % % Diagnostics tools: to make sure the gradient and Hessian are
    % % correct during the prototyping stage.
    % checkgradient(problem); pause;
    % checkhessian(problem); pause;
    
    % % To investigate the effect of the rotational invariance when using
    % % the oblique or the elliptope geometry, or to study the saddle point
    % % issue mentioned above, it is sometimes interesting to look at the
    % % spectrum of the Hessian. For large dimensions, this is slow!
    % stairs(sort(hessianspectrum(problem, Y)));
    % drawnow; pause;
    
    
    % % When facing a saddle point issue as described in the master
    % % function, and when no sure mechanism exists to find an escape
    % % direction, it may be helpful to set useRand to true and raise
    % % miniter to more than 1, when using trustregions. This will tell the
    % % solver to not stop before at least miniter iterations were
    % % accomplished (thus disregarding the zero gradient at the saddle
    % % point) and to use random search directions to kick start the inner
    % % solve (tCG) step. It is not as efficient as finding a sure escape
    % % direction, but sometimes it's the best we have.
    % options.useRand = true;
    % options.miniter = 5;
    
    options.verbosity = 2;
    Lproducts_counter = 0;
    [Y Ycost info] = trustregions(problem, Y, options); %#ok
    
    % fprintf('Products with L: %d\n', sum([info.Lproducts]));

end
