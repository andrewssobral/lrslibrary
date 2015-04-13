function candoit = canGetGradient(problem)
% Checks whether the gradient can be computed for a problem structure.
% 
% function candoit = canGetGradient(problem)
%
% Returns true if the gradient of the cost function can be computed given
% the problem description, false otherwise.
%
% See also: canGetCost canGetDirectionalDerivative canGetHessian

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 

    candoit = isfield(problem, 'grad') || isfield(problem, 'costgrad') || ...
              canGetEuclideanGradient(problem);
    
end
