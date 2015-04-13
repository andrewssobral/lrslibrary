function [stepsize, newx, storedb, lsmem, lsstats] = ...
           linesearch_hint(problem, x, d, f0, df0, options, storedb, lsmem)
% Armijo line-search based on the line-search hint in the problem structure.
%
% function [stepsize, newx, storedb, lsmem, lsstats] = 
%          linesearch_hint(problem, x, d, f0, df0, options, storedb, lsmem)
%
% Base line-search algorithm for descent methods, based on a simple
% backtracking method. The search direction provided has to be a descent
% direction, as indicated by a negative df0 = directional derivative of f
% at x along d.
%
% The algorithm obtains an initial step size candidate from the problem
% structure, typically through the problem.linesearch function. If that
% step does not fulfill the Armijo sufficient decrease criterion, that step
% size is reduced geometrically until a satisfactory step size is obtained
% or until a failure criterion triggers.
% 
% Below, the step will be constructed as alpha*d, and the step size is the
% norm of that vector, thus: stepsize = alpha*norm_d. The step is executed
% by retracting the vector alpha*d from the current point x, giving newx.
%
% Inputs
%
%  problem : structure holding the description of the optimization problem
%  x : current point on the manifold problem.M
%  d : tangent vector at x (descent direction)
%  f0 : cost value at x
%  df0 : directional derivative at x along d
%  options : options structure (see in code for usage)
%  storedb : store database structure for caching purposes
%  lsmem : not used
%
% Outputs
%
%  stepsize : norm of the vector retracted to reach newx from x.
%  newx : next iterate suggested by the line-search algorithm, such that
%         the retraction at x of the vector alpha*d reaches newx.
%  storedb : the (possibly updated) store database structure.
%  lsmem : not used.
%  lsstats : statistics about the line-search procedure (stepsize, number
%            of trials etc).
%
% See also: steepestdescent conjugategradients linesearch

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, July 17, 2014.
% Contributors: 
% Change log: 
%


    % Backtracking default parameters. These can be overwritten in the
    % options structure which is passed to the solver.
    default_options.ls_contraction_factor = .5;
    default_options.ls_suff_decr = 1e-4;
    default_options.ls_max_steps = 25;
    
    options = mergeOptions(default_options, options);
    
    contraction_factor = options.ls_contraction_factor;
    suff_decr = options.ls_suff_decr;
    max_ls_steps = options.ls_max_steps;
    
    % Obtain an initial guess at alpha from the problem structure. It is
    % assumed that the present line-search is only called when the problem
    % structure provides enough information for the call here to work.
    [alpha, storedb] = getLinesearch(problem, x, d, storedb);
    
    % Make the chosen step and compute the cost there.
    newx = problem.M.retr(x, d, alpha);
    [newf, storedb] = getCost(problem, newx, storedb);
    cost_evaluations = 1;
    
    % Backtrack while the Armijo criterion is not satisfied
    while newf > f0 + suff_decr*alpha*df0
        
        % Reduce the step size,
        alpha = contraction_factor * alpha;
        
        % and look closer down the line
        newx = problem.M.retr(x, d, alpha);
        [newf, storedb] = getCost(problem, newx, storedb);
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
    norm_d = problem.M.norm(x, d);
    stepsize = alpha * norm_d;
    
    % Save some statistics also, for possible analysis.
    lsstats.costevals = cost_evaluations;
    lsstats.stepsize = stepsize;
    lsstats.alpha = alpha;
    
end
