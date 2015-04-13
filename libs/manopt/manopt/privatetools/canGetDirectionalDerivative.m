function candoit = canGetDirectionalDerivative(problem)
% Checks whether dir. derivatives can be computed for a problem structure.
% 
% function candoit = canGetDirectionalDerivative(problem)
%
% Returns true if the directional derivatives of the cost function can be
% computed given the problem description, false otherwise.
%
% See also: canGetCost canGetGradient canGetHessian

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 

    candoit = isfield(problem, 'diff') || canGetGradient(problem);
    
end
