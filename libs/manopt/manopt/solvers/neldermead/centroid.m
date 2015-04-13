function y = centroid(M, x)
% Attempts the computation of a centroid of a set of points on amanifold.
% 
% function y = centroid(M, x)
%
% M is a structure representing a manifold. x is a cell of points on that
% manifold.

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 


    % For now, just apply a few steps of gradient descent for Karcher means
    
    n = numel(x);
    
    problem.M = M;
    
    problem.cost = @cost;
    function val = cost(y)
        val = 0;
        for i = 1 : n
            val = val + M.dist(y, x{i})^2;
        end
        val = val/2;
    end

    problem.grad = @grad;
    function g = grad(y)
        g = M.zerovec(y);
        for i = 1 : n
            g = M.lincomb(y, 1, g, -1, M.log(y, x{i}));
        end
    end

    % This line can be uncommented to check that the gradient is indeed
    % correct. This should always be the case if the dist and the log
    % functions in the manifold are correct.
    % checkgradient(problem);
    
    query = warning('query', 'manopt:getHessian:approx');
    warning('off', 'manopt:getHessian:approx')
    options.verbosity = 0;
    options.maxiter = 15;
    y = trustregions(problem, x{randi(n)}, options);
    warning(query.state, 'manopt:getHessian:approx')

end
