function [x cutvalue cutvalue_upperbound Y] = maxcut_octave(L, r)
% Algorithm to (try to) compute a maximum cut of a graph, via SDP approach.
% 
% function x = maxcut_octave(L)
% function [x cutvalue cutvalue_upperbound Y] = maxcut_octave(L, r)
%
% See examples/maxcut.m for help about the math behind this example. This
% file is here to illustrate how to use Manopt within Octave.
%
% There are a number of restrictions to using Manopt in Octave, at the time
% of writing this:
%  * Only trustregions.m works as a solver yet.
%  * Only elliptopefactory.m works as a manifold factory yet.
%  * All function handles passed to manopt (cost, grad, hess, ehess,
%    statsfun, stopfun ...) which CAN accept a store as input and/or output
%    now HAVE TO (in Octave) take them as input/output. Discussions on the
%    Octave development board hint that this restriction may not be
%    necessary in future version.
%  * You cannot define those functions as nested functions. Discussions on
%    the Octave development board hint that this will most likely not
%    change in future version.
%
% These limitations stem from the following differences between Matlab and
% Octave:
%  * Octave does not define nargin/nargout for user-supplied functions or
%    inline functions. This will likely change.
%  * Octave has no nested functions support. This will likely not change.
% Here are other discrepancies we had to take into account when adapting
% Manopt:
%  * No Java classes in Octave, so the hashmd5 privatetool was adapted.
%  * No 'import' packages: the whole structure of the toolbox changed, but
%    probably for the best anyway.
%  * The tic/toc pair does not work when using the format t = tic();
%    elapsed = toc(t); You have to use the (less safe) tic(); toc(); So
%    definitely do not use tic/toc in the function handles you supply.
%  * try/catch blocks do not give the catch an exception object.
%  * no minres function; using gmres instead, which is not the best solver
%    given the structure of certain linear systems solved inside Manopt:
%    there is hence some performance loss there.
%
% See also: maxcut

% This file is part of Manopt and is copyrighted. See the license file.
%
% Main author: Nicolas Boumal, Aug. 22, 2013
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
            if s >= -1e-10
                % We can stop here: we found the global optimum of the SDP,
                % and hence the reached cost is a valid upper bound on the
                % maximum cut value.
                cutvalue_upperbound = max(-[info.cost]);
                break;
            end
            
            % This is our escape direction.
            Z = manifold.proj(Y0, [zeros(n, rr-1) v]);
            
            % % These instructions can be uncommented to see what the cost
            % % function looks like at a saddle point.
            % plotprofile(problem, Y0, Z, linspace(-1, 1, 101));
            % drawnow; pause;
            
            % Now make a step in the Z direction to escape from the saddle.
            % It is not obvious that it is ok to do a unit step ... perhaps
            % need to be cautious here with the stepsize. It's not too
            % critical though: the important point is to leave the saddle
            % point. But it's nice to guarantee monotone decrease of the
            % cost, and we can't do that with a constant step (at least,
            % not without a proper argument to back it up).
            stepsize = 1.0;
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
    
    % % Unfortunately, you cannot code things this way in Octave, because
    % you have to accept the store as input AND return it as second output.
    % problem.cost = @(Y)  -trace(Y'*L*Y)/4;
    % problem.egrad = @(Y) -(L*Y)/2;
    % problem.ehess = @(Y, U) -(L*U)/2;
    
    % Instead of the prototyping version, the functions below describe the
    % cost, gradient and Hessian using the caching system (the store
    % structure). This alows to execute exactly the required number of
    % multiplications with the matrix L.

    problem.cost = @(Y, store) cost(L, Y, store);

    problem.grad = @(Y, store) grad(manifold, L, Y, store);

    problem.hess = @(Y, U, store) hess(manifold, L, Y, U, store);    

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
    % profile clear; profile on;
    [Y Ycost info] = trustregions(problem, Y, options); %#ok
    % profile off; profile report;

end


function store = prepare(L, Y, store)
    if ~isfield(store, 'LY')
        store.LY = L*Y;
    end
end

function [f store] = cost(L, Y, store)
    store = prepare(L, Y, store);
    LY = store.LY;
    f = -(Y(:)'*LY(:))/4; % = -trace(Y'*LY)/4;
end

function [g store] = grad(manifold, L, Y, store)
    store = prepare(L, Y, store);
    LY = store.LY;
    g = manifold.egrad2rgrad(Y, -LY/2);
end

function [h store] = hess(manifold, L, Y, U, store)
    store = prepare(L, Y, store);
    LY = store.LY;
    LU = L*U;
    h = manifold.ehess2rhess(Y, -LY/2, -LU/2, U);
end
