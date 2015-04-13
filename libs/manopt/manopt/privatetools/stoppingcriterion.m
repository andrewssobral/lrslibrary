function [stop reason] = stoppingcriterion(problem, x, options, info, last)
% Checks for standard stopping criteria, as a helper to solvers.
%
% function [stop reason] = stoppingcriterion(problem, x, options, info, last)
%
% Executes standard stopping criterion checks, based on what is defined in
% the info(last) stats structure and in the options structure.
%
% The returned number 'stop' is 0 if none of the stopping criteria
% triggered, and a (strictly) positive integer otherwise. The integer
% identifies which criterion triggered:
%  0 : Nothing triggered;
%  1 : Cost tolerance reached;
%  2 : Gradient norm tolerance reached;
%  3 : Max time exceeded;
%  4 : Max iteration count reached;
%  5 : Maximum number of cost evaluations reached;
%  6 : User defined stopfun criterion triggered.
%
% The output 'reason' is a string describing the triggered event.

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 


    stop = 0;
    reason = '';
    
    stats = info(last);

    % Target cost attained
    if isfield(stats, 'cost') && isfield(options, 'tolcost') && ...
       stats.cost <= options.tolcost
        reason = 'Cost tolerance reached. See options.tolcost.';
        stop = 1;
        return;
    end

    % Target gradient norm attained
    if isfield(stats, 'gradnorm') && isfield(options, 'tolgradnorm') && ...
       stats.gradnorm < options.tolgradnorm
        reason = 'Gradient norm tolerance reached. See options.tolgradnorm.';
        stop = 2;
        return;
    end

    % Alloted time exceeded
    if isfield(stats, 'time') && isfield(options, 'maxtime') && ...
       stats.time >= options.maxtime
        reason = 'Max time exceeded. See options.maxtime.';
        stop = 3;
        return;
    end

    % Alloted iteration count exceeded
    if isfield(stats, 'iter') && isfield(options, 'maxiter') && ...
       stats.iter >= options.maxiter
        reason = 'Max iteration count reached. See options.maxiter.';
        stop = 4;
        return;
    end
    
    % Alloted function evaluation count exceeded
    if isfield(stats, 'costevals') && isfield(options, 'maxcostevals') && ...
       stats.costevals >= options.maxcostevals
        reason = 'Maximum number of cost evaluations reached. See options.maxcostevals.';
        stop = 5;
    end

    % Check whether the possibly user defined stopping criterion
    % triggers or not.
    if isfield(options, 'stopfun')
        userstop = options.stopfun(problem, x, info, last);
        if userstop
            reason = 'User defined stopfun criterion triggered. See options.stopfun.';
            stop = 6;
            return;
        end
    end

end
