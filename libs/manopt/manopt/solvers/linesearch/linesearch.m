function [stepsize, newx, storedb, lsmem, lsstats] = ...
                linesearch(problem, x, d, f0, df0, options, storedb, lsmem)
% Standard line-search algorithm (step size selection) for descent methods.
%
% function [stepsize, newx, storedb, lsmem, lsstats] = 
%               linesearch(problem, x, d, f0, df0, options, storedb, lsmem)
%
% Base line-search algorithm for descent methods, based on a simple
% backtracking method. The search direction provided has to be a descent
% direction, as indicated by a negative df0 = directional derivative of f
% at x along d.
%
% The algorithm is invariant under positive scaling of the cost function
% and under offsetting, that is: if the cost function f is replaced by
% 8*f+3 for example, the returned step size will be the same. Furthermore,
% the returned step size is independent of the norm of the search direction
% vector d: only the direction of d is important.
% 
% Below, the step will be constructed as alpha*d, and the step size is the
% norm of that vector, thus: stepsize = alpha*norm_d. The step is executed
% by retracting the vector alpha*d from the current point x, giving newx.
%
% Inputs
%
%  problem : structure holding the description of the optimization problem
%  x : current point on the manifold problem.M
%  d : tangent vector at x (descent direction) -- its norm is irrelevant
%  f0 : cost value at x
%  df0 : directional derivative at x along d
%  options : options structure (see in code for usage)
%  storedb : store database structure for caching purposes
%  lsmem : a special memory holder to give the linesearch the opportunity
%          to "remember" what happened in the previous calls. See below.
%
% Outputs
%
%  stepsize : norm of the vector retracted to reach newx from x.
%  newx : next iterate suggested by the line-search algorithm, such that
%         the retraction at x of the vector alpha*d reaches newx.
%  storedb : the (possibly updated) store database structure.
%  lsmem : the (possibly updated) lsmem memory holder.
%  lsstats : statistics about the line-search procedure (stepsize, number
%            of trials etc).
% 
% About lsmem : It can be anything, and will typically be a matrix or a
%               structure. When first calling the line-search, it should be
%               passed as the empty matrix []. For subsequent calls
%               (pertaining to the same solver call), the previously
%               returned lsmem should be passed as input. This memory
%               holder gives the line-search a chance to exploit knowledge
%               of previous decisions to make new decisions.
%
% See also: steepestdescent conjugategradients linesearch_adaptive

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 
%     Sept. 13, 2013 (NB) : User control over the parameters of the
%                           linesearch via the options
%                           ls_contraction_factor, ls_optimism,
%                           ls_suff_decr and ls_max_steps. See in code for
%                           the effect of those.
%     Sept. 13, 2013 (NB) : The automatic direction reversal feature was
%                           removed (it triggered when df0 > 0). Direction
%                           reversal is a decision that needs to be made by
%                           the solver, so it can know about it.
%     Sept. 13, 2013 (NB) : The linesearch is now invariant under rescaling
%                           of the cost function f. That is, if f is
%                           replaced by 8*f (and hence the directional
%                           derivatives of f are scaled accordingly), the
%                           stepsizes computed will not change.
%     Nov. 7, 2013 (NB) :   The linesearch is now invariant under rescaling
%                           of the search direction d. The meaning of
%                           stepsize is also more clear in the comments.
%                           Added a parameter ls_initial_stepsize to give
%                           users control over the first step size trial.


    % Backtracking default parameters. These can be overwritten in the
    % options structure which is passed to the solver.
    default_options.ls_contraction_factor = .5;
    default_options.ls_optimism = 1/.5;
    default_options.ls_suff_decr = 1e-4;
    default_options.ls_max_steps = 25;
    default_options.ls_initial_stepsize = 1;
    
    options = mergeOptions(default_options, options);
    
    contraction_factor = options.ls_contraction_factor;
    optimism = options.ls_optimism;
    suff_decr = options.ls_suff_decr;
    max_ls_steps = options.ls_max_steps;
    initial_stepsize = options.ls_initial_stepsize;
    
    % Compute the norm of the search direction.
    % This is useful to make the linesearch algorithm invariant under the
    % scaling of d. The rationale is that the important information is the
    % search direction, not the size of that vector. The question of how
    % far we should go is precisely what the linesearch algorithm is
    % supposed to answer: the calling algorithm should not need to care.
    norm_d = problem.M.norm(x, d);
    
    % At first, we have no idea of what the step size should be.
    alpha = NaN;
    
    % If we know about what happened at the previous step, we can leverage
    % that to compute an initial guess for the step size, as inspired from
    % Nocedal&Wright, p59.
    if isstruct(lsmem) && isfield(lsmem, 'f0')
        % Pick initial step size based on where we were last time,
        alpha = 2*(f0-lsmem.f0)/df0;
        % and go look a little further (or less far), just in case.
        alpha = optimism*alpha;
    end
    
    % If we have no information about the previous iteration (maybe this is
    % the first one?) or if the above formula gave a too small step size
    % (perhaps it is even negative), then fall back to a user supplied
    % suggestion for the first step size (the "a priori").
    % At any rate, the choice should be invariant under rescaling of the
    % cost function f and of the search direction d, and it should be
    % bounded away from zero for convergence guarantees. We must allow it
    % to be close to zero though, for fine convergence.
    if isnan(alpha) || alpha*norm_d <= eps
        alpha = initial_stepsize/norm_d;
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

    % Save the situtation faced now so that at the next iteration, we will
    % know something about the previous decision.
    lsmem.f0 = f0;
    lsmem.df0 = df0;
    lsmem.stepsize = stepsize;
    
    % Save some statistics also, for possible analysis.
    lsstats.costevals = cost_evaluations;
    lsstats.stepsize = stepsize;
    lsstats.alpha = alpha;
    
end
