function checkdiff(problem, x, d)
% Checks the consistency of the cost function and directional derivatives.
%
% function checkdiff(problem)
% function checkdiff(problem, x)
% function checkdiff(problem, x, d)
%
% checkdiff performs a numerical test to check that the directional
% derivatives defined in the problem structure agree up to first order with
% the cost function at some point x, along some direction d. The test is
% based on a truncated Taylor series (see online Manopt documentation).
%
% Both x and d are optional and will be sampled at random if omitted.
%
% See also: checkgradient checkhessian

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 

        
    % Verify that the problem description is sufficient.
    if ~canGetCost(problem)
        error('It seems no cost was provided.');  
    end
    if ~canGetDirectionalDerivative(problem)
        error('It seems no directional derivatives were provided.');    
    end
    
    dbstore = struct();
        
    x_isprovided = exist('x', 'var') && ~isempty(x);
    d_isprovided = exist('d', 'var') && ~isempty(d);
    
    if ~x_isprovided && d_isprovided
        error('If d is provided, x must be too, since d is tangent at x.');
    end
    
    % If x and / or d are not specified, pick them at random.
    if ~x_isprovided
        x = problem.M.rand();
    end
    if ~d_isprovided
        d = problem.M.randvec(x);
    end

    % Compute the value f0 at f and directional derivative at x along d.
    f0 = getCost(problem, x, dbstore);
    df0 = getDirectionalDerivative(problem, x, d, dbstore);
    
    % Compute the value of f at points on the geodesic (or approximation of
    % it) originating from x, along direction d, for stepsizes in a large
    % range given by h.
    h = logspace(-8, 0, 51);
    value = zeros(size(h));
    for i = 1 : length(h)
        y = problem.M.exp(x, d, h(i));
        value(i) = getCost(problem, y, dbstore);
    end
    
    % Compute the linear approximation of the cost function using f0 and
    % df0 at the same points.
    model = polyval([df0 f0], h);
    
    % Compute the approximation error
    err = abs(model - value);
    
    % And plot it.
    loglog(h, err);
    title(sprintf(['Directional derivative check.\nThe slope of the '...
                   'continuous line should match that of the dashed '...
                   '(reference) line\nover at least a few orders of '...
                   'magnitude for h.']));
    xlabel('h');
    ylabel('Approximation error');
    
    line('xdata', [1e-8 1e0], 'ydata', [1e-8 1e8], ...
         'color', 'k', 'LineStyle', '--', ...
         'YLimInclude', 'off', 'XLimInclude', 'off');
    
    
    % In a numerically reasonable neighborhood, the error should decrease
    % as the square of the stepsize, i.e., in loglog scale, the error
    % should have a slope of 2.
    window_len = 10;
    [range poly] = identify_linear_piece(log10(h), log10(err), window_len);
    hold on;
        loglog(h(range), 10.^polyval(poly, log10(h(range))), ...
               'r-', 'LineWidth', 3);
    hold off;
    
    fprintf('The slope should be 2. It appears to be: %g.\n', poly(1));
    fprintf(['If it is far from 2, then directional derivatives ' ...
             'might be erroneous.\n']);
    
end
