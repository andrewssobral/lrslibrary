function cost = plotprofile(problem, x, d, t)
% Plot the cost function along a geodesic or a retraction path.
%
% function plotprofile(problem, x, d, t)
% function costs = plotprofile(problem, x, d, t)
%
% Plot profile evaluates the cost function along a geodesic gamma(t) such
% that gamma(0) = x and the derivative of gamma at 0 is the direction d.
% The input t is a vector specifying for which values of t we must evaluate
% f(gamma(t)) (it may include negative values).
%
% If the function is called with an output, the plot is not drawn and the
% values of the cost are returned for the instants t.

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Jan. 9, 2013.
% Contributors: 
% Change log: 

    % Verify that the problem description is sufficient.
    if ~canGetCost(problem)
        error('It seems no cost was provided.');  
    end
    
    linesearch_fun = @(t) getCost(problem, problem.M.exp(x, d, t), struct());
    
    cost = zeros(size(t));
    for i = 1 : numel(t)
        cost(i) = linesearch_fun(t(i));
    end
    
    if nargout == 0
        plot(t, cost);
        xlabel('t');
        ylabel('f(Exp_x(t*d))');
    end
    
end
