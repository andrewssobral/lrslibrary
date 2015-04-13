function checkhessian(problem, x, d)
% Checks the consistency of the cost function and the Hessian.
%
% function checkhessian(problem)
% function checkhessian(problem, x)
% function checkhessian(problem, x, d)
%
% checkhessian performs a numerical test to check that the directional
% derivatives and Hessian defined in the problem structure agree up to
% second order with the cost function at some point x, along some direction
% d. The test is based on a truncated Taylor series (see online Manopt
% documentation).
% 
% It is also tested that the Hessian along some direction is indeed a
% tangent vector and that the Hessian operator is symmetric w.r.t. the
% Riemannian metric.
% 
% Both x and d are optional and will be sampled at random if omitted.
%
% See also: checkdiff checkgradient

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 

        
    % Verify that the problem description is sufficient.
    if ~canGetCost(problem)
        error('It seems no cost was provided.');  
    end
    if ~canGetGradient(problem)
        error('It seems no gradient provided.');    
    end
    if ~canGetHessian(problem)
        error('It seems no Hessian was provided.');    
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
    
    %% Check that the directional derivative and the Hessian at x along d
    %% yield a second order model of the cost function.
    
    % Compute the value f0 at f, directional derivative df0 at x along d,
    % and Hessian along [d, d].
    f0 = getCost(problem, x, dbstore);
    df0 = getDirectionalDerivative(problem, x, d, dbstore);
    d2f0 = problem.M.inner(x, d, getHessian(problem, x, d, dbstore));
    
    % Compute the value of f at points on the geodesic (or approximation of
    % it) originating from x, along direction d, for stepsizes in a large
    % range given by h.
    h = logspace(-8, 0, 51);
    value = zeros(size(h));
    for i = 1 : length(h)
        y = problem.M.exp(x, d, h(i));
        value(i) = getCost(problem, y, dbstore);
    end
    
    % Compute the quadratic approximation of the cost function using f0,
    % df0 and d2f0 at the same points.
    model = polyval([.5*d2f0 df0 f0], h);
    
    % Compute the approximation error
    err = abs(model - value);
    
    % And plot it.
    loglog(h, err);
    title(sprintf(['Hessian check.\nThe slope of the continuous line ' ...
                   'should match that of the dashed (reference) line\n' ...
                   'over at least a few orders of magnitude for h.']));
    xlabel('h');
    ylabel('Approximation error');
    
    line('xdata', [1e-8 1e0], 'ydata', [1e-16 1e8], ...
         'color', 'k', 'LineStyle', '--', ...
         'YLimInclude', 'off', 'XLimInclude', 'off');
    
    % In a numerically reasonable neighborhood, the error should decrease
    % as the cube of the stepsize, i.e., in loglog scale, the error
    % should have a slope of 3.
    window_len = 10;
    [range poly] = identify_linear_piece(log10(h), log10(err), window_len);
    hold on;
        loglog(h(range), 10.^polyval(poly, log10(h(range))), ...
               'r-', 'LineWidth', 3);
    hold off;
    
    fprintf('The slope should be 3. It appears to be: %g.\n', poly(1));
    fprintf(['If it is far from 3, then directional derivatives or ' ...
             'the Hessian might be erroneous.\n']);

    
    %% Check that the Hessian at x along direction d is a tangent vector.
    if isfield(problem.M, 'tangent')
        hess = getHessian(problem, x, d, dbstore);
        phess = problem.M.tangent(x, hess);
        residual = problem.M.lincomb(x, 1, hess, -1, phess);
        err = problem.M.norm(x, residual);
        fprintf('The residual should be zero, or very close. ');
        fprintf('Residual: %g.\n', err);
        fprintf(['If it is far from 0, then the Hessian is not in the ' ...
                 'tangent plane.\n']);
    else
        fprintf(['Unfortunately, Manopt was un able to verify that the '...
                 'Hessian is indeed a tangent vector.\nPlease verify ' ...
                 'this manually.']);
    end    
    
    %% Check that the Hessian at x is symmetric.
    d1 = problem.M.randvec(x);
    d2 = problem.M.randvec(x);
    h1 = getHessian(problem, x, d1, dbstore);
    h2 = getHessian(problem, x, d2, dbstore);
    v1 = problem.M.inner(x, d1, h2);
    v2 = problem.M.inner(x, h1, d2);
    value = v1-v2;
    fprintf(['<d1, H[d2]> - <H[d1], d2> should be zero, or very close.' ...
             '\n\tValue: %g - %g = %g.\n'], v1, v2, value);
    fprintf('If it is far from 0, then the Hessian is not symmetric.\n');
    
end
