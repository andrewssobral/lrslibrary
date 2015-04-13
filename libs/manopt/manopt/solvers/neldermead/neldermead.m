function [x, cost, info, options] = neldermead(problem, x, options)
% Nelder Mead optimization algorithm for derivative-free minimization.
%
% function [x, cost, info, options] = neldermead(problem)
% function [x, cost, info, options] = neldermead(problem, x0)
% function [x, cost, info, options] = neldermead(problem, x0, options)
% function [x, cost, info, options] = neldermead(problem, [], options)
%
% Apply the Nelder-Mead minimization algorithm to the problem defined in
% the problem structure, starting with the population x0 if it is provided
% (otherwise, a random population on the manifold is generated). A
% population is a cell containing points on the manifold. The number of
% elements in the cell must be dim+1, where dim is the dimension of the
% manifold: problem.M.dim().
%
% To specify options whilst not specifying an initial guess, give x0 as []
% (the empty matrix).
%
% This algorithm is a plain adaptation of the Euclidean Nelder-Mead method
% to the Riemannian setting. It comes with no convergence guarantees and
% there is room for improvement. In particular, we compute centroids as
% Karcher means, which seems overly expensive: cheaper forms of
% average-like quantities might work better.
% This solver is useful nonetheless for problems for which no derivatives
% are available, and it may constitute a starting point for the development
% of other Riemannian derivative-free methods.
%
% None of the options are mandatory. See in code for details.
%
% Based on http://www.optimization-online.org/DB_FILE/2007/08/1742.pdf.
%
% See also: manopt/solvers/pso/pso

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 

    
    % Verify that the problem description is sufficient for the solver.
    if ~canGetCost(problem)
        warning('manopt:getCost', ...
                'No cost provided. The algorithm will likely abort.');  
    end
    
    % Dimension of the manifold
    dim = problem.M.dim();

    % Set local defaults here
    localdefaults.storedepth = 0;                     % no need for caching
    localdefaults.maxcostevals = max(1000, 2*dim);
    localdefaults.maxiter = max(2000, 4*dim);
    
    localdefaults.reflection = 1;
    localdefaults.expansion = 2;
    localdefaults.contraction = .5;
    % forced to .5 to enable using pairmean functions in manifolds.
    % localdefaults.shrinkage = .5;
    
    % Merge global and local defaults, then merge w/ user options, if any.
    localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    options = mergeOptions(localdefaults, options);
    
    % Create a store database.
    storedb = struct();
    
    % Start timing for initialization.
    timetic = tic();
    
    % If no initial simplex x is given by the user, generate one at random.
    if ~exist('x', 'var') || isempty(x)
        x = cell(dim+1, 1);
        for i = 1 : dim+1
            x{i} = problem.M.rand();
        end
    end
    
    % Compute objective-related quantities for x, and setup a
    % function evaluations counter.
    costs = zeros(dim+1, 1);
    for i = 1 : dim+1
        [costs(i) storedb] = getCost(problem, x{i}, storedb);
    end
    costevals = dim+1;
    
    % Sort simplex points by cost.
    [costs order] = sort(costs);
    x = x(order);
    
    % Iteration counter (at any point, iter is the number of fully executed
    % iterations so far).
    iter = 0;
    
    % Save stats in a struct array info, and preallocate
    % (see http://people.csail.mit.edu/jskelly/blog/?x=entry:entry091030-033941)
    % savestats will be called twice for the initial iterate (number 0),
    % which is unfortunate, but not problematic.
    stats = savestats();
    info(1) = stats;
    info(min(10000, options.maxiter+1)).iter = [];
    
    % Start iterating until stopping criterion triggers.
    while true
        
        % Make sure we don't use to much memory for the store database.
        storedb = purgeStoredb(storedb, options.storedepth);
        
        stats = savestats();
        info(iter+1) = stats; %#ok<AGROW>
        iter = iter + 1;
        
        % Start timing this iteration.
        timetic = tic();
        
        % Sort simplex points by cost.
        [costs order] = sort(costs);
        x = x(order);

        % Log / display iteration information here.
        if options.verbosity >= 2
            fprintf('Cost evals: %7d\tBest cost: %+.4e\t', ...
                    costevals, costs(1));
        end
        
        % Run standard stopping criterion checks.
        [stop reason] = stoppingcriterion(problem, x, options, info, iter);
    
        if stop
            if options.verbosity >= 1
                fprintf([reason '\n']);
            end
            break;
        end
        
        % Compute a centroid for the dim best points.
        xbar = centroid(problem.M, x(1:end-1));
        
        % Compute the direction for moving along the axis xbar - worst x.
        vec = problem.M.log(xbar, x{end});
        
        % Reflection step
        xr = problem.M.exp(xbar, vec, -options.reflection);
        [costr storedb] = getCost(problem, xr, storedb);
        costevals = costevals + 1;
        
        % If the reflected point is honorable, drop the worst point,
        % replace it by the reflected point and start new iteration.
        if costr >= costs(1) && costr < costs(end-1)
            fprintf('Reflection\n');
            costs(end) = costr;
            x{end} = xr;
            continue;
        end
        
        % If the reflected point is better than our best point, expand.
        if costr < costs(1)
            xe = problem.M.exp(xbar, vec, -options.expansion);
            [coste storedb] = getCost(problem, xe, storedb);
            costevals = costevals + 1;
            if coste < costr
                fprintf('Expansion\n');
                costs(end) = coste;
                x{end} = xe;
                continue;
            else
                fprintf('Reflection (failed expansion)\n');
                costs(end) = costr;
                x{end} = xr;
                continue;
            end
        end
        
        % If the reflected point is worse than the second to worst point,
		% contract.
        if costr >= costs(end-1)
            if costr < costs(end)
                % do an outside contraction
                xoc = problem.M.exp(xbar, vec, -options.contraction);
                [costoc storedb] = getCost(problem, xoc, storedb);
                costevals = costevals + 1;
                if costoc <= costr
                    fprintf('Outside contraction\n');
                    costs(end) = costoc;
                    x{end} = xoc;
                    continue;
                end
            else
                % do an inside contraction
                xic = problem.M.exp(xbar, vec, options.contraction);
                [costic storedb] = getCost(problem, xic, storedb);
                costevals = costevals + 1;
                if costic <= costs(end)
                    fprintf('Inside contraction\n');
                    costs(end) = costic;
                    x{end} = xic;
                    continue;
                end
            end
        end
        
        % If we get here, shrink the simplex around x{1}.
        fprintf('Shrinkage\n');
        for i = 2 : dim+1
            x{i} = problem.M.pairmean(x{1}, x{i});
            [costs(i) storedb] = getCost(problem, x{i}, storedb);
        end
        costevals = costevals + dim;
        
    end
    
    
    info = info(1:iter);
    cost = costs(1);
    x = x{1};
    
    
    
    
    % Routine in charge of collecting the current iteration stats.
    function stats = savestats()
        stats.iter = iter;
        stats.cost = costs(1);
        stats.costevals = costevals;
        if iter == 0
            stats.time = toc(timetic);
        else
            stats.time = info(iter).time + toc(timetic);
        end
        stats = applyStatsfun(problem, x, storedb, options, stats);
    end
    
end
