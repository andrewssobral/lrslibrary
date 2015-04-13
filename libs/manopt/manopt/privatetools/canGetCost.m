function candoit = canGetCost(problem)
% Checks whether the cost function can be computed for a problem structure.
%
% function candoit = canGetCost(problem)
%
% Returns true if the cost function can be computed given the problem
% description, false otherwise.
%
% See also: getCost canGetDirectionalDerivative canGetGradient canGetHessian

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 


    candoit = isfield(problem, 'cost') || isfield(problem, 'costgrad');
    
end
