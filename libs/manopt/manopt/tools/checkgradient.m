function checkgradient(problem, x, d)
% Checks the consistency of the cost function and the gradient.
%
% function checkgradient(problem)
% function checkgradient(problem, x)
% function checkgradient(problem, x, d)
%
% checkgradient performs a numerical test to check that the gradient
% defined in the problem structure agrees up to first order with the cost
% function at some point x, along some direction d. The test is based on a
% truncated Taylor series (see online Manopt documentation).
%
% It is also tested that the gradient is indeed a tangent vector.
% 
% Both x and d are optional and will be sampled at random if omitted.
%
% See also: checkdiff checkhessian

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

    %% Check that the gradient yields a first order model of the cost.
    
    % By removing the 'directional derivative' function, it should be so
    % (?) that the checkdiff function will use the gradient to compute
    % directional derivatives.
    if isfield(problem, 'diff')
        problem = rmfield(problem, 'diff');
    end
    checkdiff(problem, x, d);
    title(sprintf(['Gradient check.\nThe slope of the continuous line ' ...
                   'should match that of the dashed (reference) line\n' ...
                   'over at least a few orders of magnitude for h.']));
    xlabel('h');
    ylabel('Approximation error');
    
    %% Try to check that the gradient is a tangent vector.
    if isfield(problem.M, 'tangent')
        grad = getGradient(problem, x, dbstore);
        pgrad = problem.M.tangent(x, grad);
        residual = problem.M.lincomb(x, 1, grad, -1, pgrad);
        err = problem.M.norm(x, residual);
        fprintf('The residual should be 0, or very close. Residual: %g.\n', err);
        fprintf('If it is far from 0, then the gradient is not in the tangent space.\n');
    else
        fprintf(['Unfortunately, Manopt was unable to verify that the '...
                 'gradient is indeed a tangent vector.\nPlease verify ' ...
                 'this manually or implement the ''tangent'' function ' ...
                 'in your manifold structure.']);
    end

end
