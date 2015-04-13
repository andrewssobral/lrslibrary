function [stepsize newx storedb lsmem lsstats] = ...
       linesearch_adaptive(problem, x, d, f0, df0, options, storedb, lsmem)
% Adaptive line search algorithm (step size selection) for descent methods.
%
% function [stepsize newx storedb lsmem lsstats] = 
%      linesearch_adaptive(problem, x, d, f0, df0, options, storedb, lsmem)
%
% Adaptive linesearch algorithm for descent methods, based on a simple
% backtracking method. On average, this line search intends to do only one
% or two cost evaluations.
%
% Contrary to linesearch.m, this function is not invariant under rescaling
% of the search direction d. Nevertheless, it sometimes performs better.
%
% Inputs/Outputs : see help for linesearch
%
% See also: steepestdescent conjugategradients linesearch

% This file is part of Manopt: www.manopt.org.
% Original author: Bamdev Mishra, Dec. 30, 2012.
% Contributors: Nicolas Boumal
% Change log:
%
%   Sept. 13, 2013 (NB) :
%       The automatic direction reversal feature was removed (it triggered
%       when df0 > 0). Direction reversal is a decision that needs to be
%       made by the solver, so it can know about it.
%
%	Nov. 7, 2013 (NB) :
%       The whole function has been recoded to mimick more closely the new
%       version of linesearch.m. The parameters are available through the
%       options structure passed to the solver and have the same names and
%       same meaning as for the base linesearch. The information is logged
%       more reliably.


    % Backtracking default parameters. These can be overwritten in the
    % options structure which is passed to the solver.
    default_options.ls_contraction_factor = .5;
    default_options.ls_suff_decr = .5;
    default_options.ls_max_steps = 10;
    default_options.ls_initial_stepsize = 1;
    options = mergeOptions(default_options, options);
    
    contraction_factor = options.ls_contraction_factor;
    suff_decr = options.ls_suff_decr;
    max_ls_steps = options.ls_max_steps;
    initial_stepsize = options.ls_initial_stepsize;
    
    % Compute the norm of the search direction.
    norm_d = problem.M.norm(x, d);
    
    % If this is not the first iteration, then lsmem should have been
    % filled with a suggestion for the initial step.
    if isstruct(lsmem) && isfield(lsmem, 'init_alpha')
        % Pick initial step size based on where we were last time,
        alpha = lsmem.init_alpha;
    
    % Otherwise, fall back to a user supplied suggestion.
    else
        alpha = initial_stepsize / norm_d;
    end

    % Make the chosen step and compute the cost there.
    newx = problem.M.retr(x, d, alpha);
    [newf storedb] = getCost(problem, newx, storedb);
    cost_evaluations = 1;
    
    % Backtrack while the Armijo criterion is not satisfied
    while newf > f0 + suff_decr*alpha*df0
        
        % Reduce the step size,
        alpha = contraction_factor * alpha;
        
        % and look closer down the line
        newx = problem.M.retr(x, d, alpha);
        [newf storedb] = getCost(problem, newx, storedb);
        cost_evaluations = cost_evaluations + 1;
        
        % Make sure we don't run out of budget
        if cost_evaluations >= max_ls_steps
            break;
        end
        
    end
    
    % If we got here without obtaining a decrease, we reject the step.
    if newf > f0
        alpha = 0;
        newx = x;
        newf = f0; %#ok<NASGU>
    end
    
    % As seen outside this function, stepsize is the size of the vector we
    % retract to make the step from x to newx. Since the step is alpha*d:
    stepsize = alpha * norm_d;

    % Fill lsmem with a suggestion for what the next initial step size
    % trial should be. On average we intend to do only one extra cost
    % evaluation. Notice how the suggestion is not about stepsize but about
    % alpha. This is the reason why this line search is not invariant under
    % rescaling of the search direction d.
    switch cost_evaluations
        case 1
            % If things go well, push your luck.
            lsmem.init_alpha = 2 * alpha;
        case 2
            % If things go smoothly, try to keep pace.
            lsmem.init_alpha = alpha;
        otherwise
            % If you backtracked a lot, the new stepsize is probably quite
            % small: try to recover.
            lsmem.init_alpha = 2 * alpha;
    end
    
    % Save some statistics also, for possible analysis.
    lsstats.costevals = cost_evaluations;
    lsstats.stepsize = stepsize;
    lsstats.alpha = alpha;
    
end
