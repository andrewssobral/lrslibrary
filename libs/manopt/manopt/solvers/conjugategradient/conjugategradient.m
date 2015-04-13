function [x, cost, info, options] = conjugategradient(problem, x, options)
% Conjugate gradient minimization algorithm for Manopt.
%
% function [x, cost, info, options] = conjugategradient(problem)
% function [x, cost, info, options] = conjugategradient(problem, x0)
% function [x, cost, info, options] = conjugategradient(problem, x0, options)
% function [x, cost, info, options] = conjugategradient(problem, [], options)
%
% Apply the conjugate gradient minimization algorithm to the problem
% defined in the problem structure, starting at x0 if it is provided
% (otherwise, at a random point on the manifold). To specify options whilst
% not specifying an initial guess, give x0 as [] (the empty matrix).
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
%   beta : value of the beta parameter (see options.beta_type)
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
%   beta_type ('H-S')
%       Conjugate gradient beta rule used to construct the new search
%       direction, based on a linear combination of the previous search
%       direction and the new (preconditioned) gradient. Possible values
%       for this parameter are:
%           'S-D', 'steep' for beta = 0 (preconditioned steepest descent)
%           'F-R' for Fletcher-Reeves's rule
%           'P-R' for Polak-Ribiere's modified rule
%           'H-S' for Hestenes-Stiefel's modified rule
%           'H-Z' for Hager-Zhang's modified rule
%       See Hager and Zhang 2006, "A survey of nonlinear conjugate gradient
%       methods" for a description of these rules in the Euclidean case and
%       for an explanation of how to adapt them to the preconditioned case.
%       The adaption to the Riemannian case is straightforward: see in code
%       for details. Modified rules take the max between 0 and the computed
%       beta value, which provides automatic restart, except for H-Z which
%       uses a different modification.
%   orth_value (Inf)
%       Following Powell's restart strategy (Math. prog. 1977), restart CG
%       (that is, make a -preconditioned- gradient step) if two successive
%       -preconditioned- gradients are "too" parallel. See for example
%       Hager and Zhang 2006, "A survey of nonlinear conjugate gradient
%       methods", page 12. An infinite value disables this strategy. See in
%       code formula for the specific criterion used.
%   linesearch (@linesearch_adaptive or @linesearch_hint)
%       Function handle to a line search function. The options structure is
%       passed to the line search too, so you can pass it parameters. See
%       each line search's documentation for info. Another available line
%       search in manopt is @linesearch, in /manopt/linesearch/linesearch.m
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
%       the CG algorithm, a store depth of 2 should always be sufficient.
%
%
% See also: steepestdescent trustregions manopt/solvers/linesearch manopt/examples

% This file is part of Manopt: www.manopt.org.
% Original author: Bamdev Mishra, Dec. 30, 2012.
% Contributors: Nicolas Boumal
% Change log: 
%
%   March 14, 2013, NB:
%       Added preconditioner support : see Section 8 in
%       https://www.math.lsu.edu/~hozhang/papers/cgsurvey.pdf
%    
%   Sept. 13, 2013, NB:
%       Now logging beta parameter too.
%    
%	Nov. 7, 2013, NB:
%       The search direction is not normalized before it is passed to the
%       linesearch anymore. This way, it is up to the designers of the
%       linesearch to decide whether they want to use the norm of the
%       search direction in their algorithm or not. There are reasons
%       against it, but practical evidence that it may help too, so we
%       allow it. The default linesearch_adaptive used does exploit the
%       norm information. The base linesearch does not. You may select it
%       by setting options.linesearch = @linesearch;
%
%	Nov. 29, 2013, NB:
%       Documentation improved: options are now explicitly described.
%       Removed the Daniel rule for beta: it was not appropriate for
%       preconditioned CG and I could not find a proper reference for it.


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
localdefaults.storedepth = 2;
% Changed by NB : H-S has the "auto restart" property.
% See Hager-Zhang 2005/2006 survey about CG methods.
% Well, the auto restart comes from the 'max(0, ...)', not so much from the
% reason stated in Hager-Zhang I believe. P-R also has auto restart.
localdefaults.beta_type = 'H-S';
localdefaults.orth_value = Inf; % by BM as suggested in Nocedal and Wright

    
% Depending on whether the problem structure specifies a hint for
% line-search algorithms, choose a default line-search that works on
% its own (typical) or that uses the hint.
if ~canGetLinesearch(problem)
    localdefaults.linesearch = @linesearch_adaptive;
else
    localdefaults.linesearch = @linesearch_hint;
end

% Merge global and local defaults, then merge w/ user options, if any.
localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
if ~exist('options', 'var') || isempty(options)
    options = struct();
end
options = mergeOptions(localdefaults, options);

% for convenience
inner = problem.M.inner;
lincomb = problem.M.lincomb;

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
[Pgrad storedb] = getPrecon(problem, x, grad, storedb);
gradPgrad = inner(x, grad, Pgrad);

% Iteration counter (at any point, iter is the number of fully executed
% iterations so far)
iter = 0;

% Save stats in a struct array info and preallocate,
% see http://people.csail.mit.edu/jskelly/blog/?x=entry:entry091030-033941
stats = savestats();
info(1) = stats;
info(min(10000, options.maxiter+1)).iter = [];

% Initial linesearch memory
lsmem = [];


if options.verbosity >= 2
    fprintf(' iter\t                cost val\t     grad. norm\n');
end

% Compute a first descent direction (not normalized)
desc_dir = lincomb(x, -1, Pgrad);


% Start iterating until stopping criterion triggers
while true
    
    % Display iteration information
    if options.verbosity >= 2
        fprintf('%5d\t%+.16e\t%.8e\n', iter, cost, gradnorm);
    end
    
    % Start timing this iteration
    timetic = tic();
    
    % Run standard stopping criterion checks
    [stop reason] = stoppingcriterion(problem, x, options, info, iter+1);
    
    % Run specific stopping criterion check
    if ~stop && abs(stats.stepsize) < options.minstepsize
        stop = true;
        reason = 'Last stepsize smaller than minimum allowed. See options.minstepsize.';
    end
    
    if stop
        if options.verbosity >= 1
            fprintf([reason '\n']);
        end
        break;
    end
    
    
    % The line search algorithms require the directional derivative of the
    % cost at the current point x along the search direction.
    df0 = inner(x, grad, desc_dir);
        
    % If we didn't get a descent direction: restart, i.e., switch to the
    % negative gradient. Equivalent to resetting the CG direction to a
    % steepest descent step, which discards the past information.
    if df0 >= 0
        
        % Or we switch to the negative gradient direction.
        if options.verbosity >= 3
            fprintf(['Conjugate gradient info: got an ascent direction '...
                     '(df0 = %2e), reset to the (preconditioned) '...
                     'steepest descent direction.\n'], df0);
        end
        % Reset to negative gradient: this discards the CG memory.
        desc_dir = lincomb(x, -1, Pgrad);
        df0 = -gradPgrad;
        
    end
    
    
    % Execute line search
    [stepsize newx storedb lsmem lsstats] = options.linesearch(...
                 problem, x, desc_dir, cost, df0, options, storedb, lsmem);

    
    % Compute the new cost-related quantities for x
    [newcost newgrad storedb] = getCostGrad(problem, newx, storedb);
    newgradnorm = problem.M.norm(newx, newgrad);
    [Pnewgrad storedb] = getPrecon(problem, x, newgrad, storedb);
    newgradPnewgrad = inner(newx, newgrad, Pnewgrad);
    
    
    % Apply the CG scheme to compute the next search direction.
    %
    % This paper https://www.math.lsu.edu/~hozhang/papers/cgsurvey.pdf
	% by Hager and Zhang lists many known beta rules. The rules defined
    % here can be found in that paper (or are provided with additional
    % references), adapted to the Riemannian setting.
	% 
    if strcmpi(options.beta_type, 'steep') || ...
       strcmpi(options.beta_type, 'S-D')              % Gradient Descent
        
        beta = 0;
        desc_dir = lincomb(x, -1, Pnewgrad);
        
    else
        
        oldgrad = problem.M.transp(x, newx, grad);
        orth_grads = inner(newx, oldgrad, Pnewgrad)/newgradPnewgrad;
        
        % Powell's restart strategy (see page 12 of Hager and Zhang's
        % survey on conjugate gradient methods, for example)
        if abs(orth_grads) >= options.orth_value,
            beta = 0;
            desc_dir = lincomb(x, -1, Pnewgrad);
            
        else % Compute the CG modification
            
            desc_dir = problem.M.transp(x, newx, desc_dir);
            
            if strcmp(options.beta_type, 'F-R')  % Fletcher-Reeves
                beta = newgradPnewgrad / gradPgrad;
                
            elseif strcmp(options.beta_type, 'P-R')  % Polak-Ribiere+
                % vector grad(new) - transported grad(current)
                diff = lincomb(newx, 1, newgrad, -1, oldgrad);
                ip_diff = inner(newx, Pnewgrad, diff);
                beta = ip_diff/gradPgrad;
                beta = max(0, beta);
                
            elseif strcmp(options.beta_type, 'H-S')  % Hestenes-Stiefel+
                diff = lincomb(newx, 1, newgrad, -1, oldgrad);
                ip_diff = inner(newx, Pnewgrad, diff);
                beta = ip_diff / inner(newx, diff, desc_dir);
                beta = max(0, beta);

            elseif strcmp(options.beta_type, 'H-Z') % Hager-Zhang+
                diff = lincomb(newx, 1, newgrad, -1, oldgrad);
                Poldgrad = problem.M.transp(x, newx, Pgrad);
                Pdiff = lincomb(newx, 1, Pnewgrad, -1, Poldgrad);
                deno = inner(newx, diff, desc_dir);
                numo = inner(newx, diff, Pnewgrad);
                numo = numo - 2*inner(newx, diff, Pdiff)*...
                                       inner(newx, desc_dir, newgrad)/deno;
                beta = numo/deno;
                
                % Robustness (see Hager-Zhang paper mentioned above)
                desc_dir_norm = problem.M.norm(newx, desc_dir);
                eta_HZ = -1/(desc_dir_norm * min(0.01, gradnorm));
                beta = max(beta,  eta_HZ);

            else
                error(['Unknown options.beta_type. ' ...
                       'Should be steep, S-D, F-R, P-R, H-S or H-Z.']);
            end
            desc_dir = lincomb(newx, -1, Pnewgrad, beta, desc_dir);
        end
        
    end
    
    % Make sure we don't use too much memory for the store database.
    storedb = purgeStoredb(storedb, options.storedepth);
    
    % Update iterate info
    x = newx;
    cost = newcost;
    grad = newgrad;
    Pgrad = Pnewgrad;
    gradnorm = newgradnorm;
    gradPgrad = newgradPnewgrad;
    
    % iter is the number of iterations we have accomplished.
    iter = iter + 1;
    
    % Log statistics for freshly executed iteration
    stats = savestats();
    info(iter+1) = stats; %#ok<AGROW>
    
end


info = info(1:iter+1);

if options.verbosity >= 1
    fprintf('Total time is %f [s] (excludes statsfun)\n', info(end).time);
end


% Routine in charge of collecting the current iteration stats
    function stats = savestats()
        stats.iter = iter;
        stats.cost = cost;
        stats.gradnorm = gradnorm;
        if iter == 0
            stats.stepsize = nan;
            stats.time = toc(timetic);
            stats.linesearch = [];
            stats.beta = 0;
        else
            stats.stepsize = stepsize;
            stats.time = info(iter).time + toc(timetic);
            stats.linesearch = lsstats;
            stats.beta = beta;
        end
        stats = applyStatsfun(problem, x, storedb, options, stats);
    end

end


