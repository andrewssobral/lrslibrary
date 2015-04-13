function [x, cost, info, options] = steepestdescent(problem, x, options)
% Steepest descent (gradient descent) minimization algorithm for Manopt.
%
% function [x, cost, info, options] = steepestdescent(problem)
% function [x, cost, info, options] = steepestdescent(problem, x0)
% function [x, cost, info, options] = steepestdescent(problem, x0, options)
% function [x, cost, info, options] = steepestdescent(problem, [], options)
%
% Apply the steepest descent minimization algorithm to the problem defined
% in the problem structure, starting at x0 if it is provided (otherwise, at
% a random point on the manifold). To specify options whilst not specifying
% an initial guess, give x0 as [] (the empty matrix).
%
% In most of the examples bundled with the toolbox (see link below), the
% solver can be replaced by the present one if need be.
%
% The outputs x and cost are the best reached point on the manifold and its
% cost. The struct-array info contains information about the iterations:
%   iter : the iteration number (0 for the initial guess)
%   cost : cost value
%   time : elapsed time in seconds
%   gradnorm : Riemannian norm of the gradient
%   stepsize : norm of the last tangent vector retracted
%   linesearch : information logged by options.linesearch
%   And possibly additional information logged by options.statsfun.
% For example, type [info.gradnorm] to obtain a vector of the successive
% gradient norms reached.
%
% The options structure is used to overwrite the default values. All
% options have a default value and are hence optional. To force an option
% value, pass an options structure with a field options.optionname, where
% optionname is one of the following and the default value is indicated
% between parentheses:
%
%   tolgradnorm (1e-6)
%       The algorithm terminates if the norm of the gradient drops below this.
%   maxiter (1000)
%       The algorithm terminates if maxiter iterations have been executed.
%   maxtime (Inf)
%       The algorithm terminates if maxtime seconds elapsed.
%   minstepsize (1e-10)
%       The algorithm terminates if the linesearch returns a displacement
%       vector (to be retracted) smaller in norm than this value.
%   linesearch (@linesearch or @linesearch_hint)
%       Function handle to a line search function. The options structure is
%       passed to the line search too, so you can pass it parameters. See
%       each line search's documentation for info. Another available line
%       search in manopt is @linesearch_adaptive, in
%       /manopt/linesearch/linesearch_adaptive.m
%       If the problem structure includes a line search hint, then the
%       default line search used in @linesearch_hint.
%   statsfun (none)
%       Function handle to a function that will be called after each
%       iteration to provide the opportunity to log additional statistics.
%       They will be returned in the info struct. See the generic Manopt
%       documentation about solvers for further information.
%   stopfun (none)
%       Function handle to a function that will be called at each iteration
%       to provide the opportunity to specify additional stopping criteria.
%       See the generic Manopt documentation about solvers for further
%       information.
%   verbosity (3)
%       Integer number used to tune the amount of output the algorithm
%       generates during execution (mostly as text in the command window).
%       The higher, the more output. 0 means silent.
%   storedepth (2)
%       Maximum number of different points x of the manifold for which a
%       store structure will be kept in memory in the storedb. If the
%       caching features of Manopt are not used, this is irrelevant. For
%       the SD algorithm, a store depth of 2 should always be sufficient.
%
%
% See also: conjugategradient trustregions manopt/solvers/linesearch manopt/examples

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 

    
    % Verify that the problem description is sufficient for the solver.
    if ~canGetCost(problem)
        warning('manopt:getCost', ...
                'No cost provided. The algorithm will likely abort.');  
    end
    if ~canGetGradient(problem)
        warning('manopt:getGradient', ...
                'No gradient provided. The algorithm will likely abort.');    
    end

    % Set local defaults here
    localdefaults.minstepsize = 1e-10;
    localdefaults.maxiter = 1000;
    localdefaults.tolgradnorm = 1e-6;
    
    % Depending on whether the problem structure specifies a hint for
    % line-search algorithms, choose a default line-search that works on
    % its own (typical) or that uses the hint.
    if ~canGetLinesearch(problem)
        localdefaults.linesearch = @linesearch;
    else
        localdefaults.linesearch = @linesearch_hint;
    end
    
    % Merge global and local defaults, then merge w/ user options, if any.
    localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    options = mergeOptions(localdefaults, options);
    
    % Create a store database
    storedb = struct();
    
    timetic = tic();
    
    % If no initial point x is given by the user, generate one at random.
    if ~exist('x', 'var') || isempty(x)
        x = problem.M.rand();
    end
    
    % Compute objective-related quantities for x
    [cost grad storedb] = getCostGrad(problem, x, storedb);
    gradnorm = problem.M.norm(x, grad);
    
    % Iteration counter (at any point, iter is the number of fully executed
    % iterations so far)
    iter = 0;
    
    % Save stats in a struct array info, and preallocate; see:
    % http://people.csail.mit.edu/jskelly/blog/?x=entry:entry091030-033941
    stats = savestats();
    info(1) = stats;
    info(min(10000, options.maxiter+1)).iter = [];
    
    % Initial line search memory
    lsmem = [];
    
    if options.verbosity >= 2
        fprintf(' iter\t                cost val\t     grad. norm\n');
    end
    
    % Start iterating until stopping criterion triggers
    while true

        % Display iteration information
        if options.verbosity >= 2
            fprintf('%5d\t%+.16e\t%.8e\n', iter, cost, gradnorm);
        end
        
        % Start timing this iteration
        timetic = tic();
        
        % Run standard stopping criterion checks
        [stop reason] = stoppingcriterion(problem, x, options, ...
                                                             info, iter+1);
        
        % If none triggered, run specific stopping criterion check
        if ~stop && stats.stepsize < options.minstepsize
            stop = true;
            reason = 'Last stepsize smaller than minimum allowed. See options.minstepsize.';
        end
    
        if stop
            if options.verbosity >= 1
                fprintf([reason '\n']);
            end
            break;
        end

        % Pick the descent direction as minus the gradient
        desc_dir = problem.M.lincomb(x, -1, grad);
        
        % Execute the line search
        [stepsize newx storedb lsmem lsstats] = options.linesearch( ...
                      problem, x, desc_dir, cost, -gradnorm^2, ...
                      options, storedb, lsmem);
        
        % Compute the new cost-related quantities for x
        [newcost newgrad storedb] = getCostGrad(problem, newx, storedb);
        newgradnorm = problem.M.norm(newx, newgrad);
        
        % Make sure we don't use too much memory for the store database
        storedb = purgeStoredb(storedb, options.storedepth);
        
        % Update iterate info
        x = newx;
        cost = newcost;
        grad = newgrad;
        gradnorm = newgradnorm;
        
        % iter is the number of iterations we have accomplished.
        iter = iter + 1;
        
        % Log statistics for freshly executed iteration
        stats = savestats();
        info(iter+1) = stats; %#ok<AGROW>
        
    end
    
    
    info = info(1:iter+1);

    if options.verbosity >= 1
        fprintf('Total time is %f [s] (excludes statsfun)\n', ...
                info(end).time);
    end
    
    
    
    % Routine in charge of collecting the current iteration stats
    function stats = savestats()
        stats.iter = iter;
        stats.cost = cost;
        stats.gradnorm = gradnorm;
        if iter == 0
            stats.stepsize = NaN;
            stats.time = toc(timetic);
            stats.linesearch = [];
        else
            stats.stepsize = stepsize;
            stats.time = info(iter).time + toc(timetic);
            stats.linesearch = lsstats;
        end
        stats = applyStatsfun(problem, x, storedb, options, stats);
    end
    
end
