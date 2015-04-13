function candoit = canGetPrecon(problem)
% Checks whether a preconditioner was specified in the problem description.
%
% function candoit = canGetPrecon(problem)
%
% Returns true if a preconditioner was specified, false otherwise. Notice
% that even if this function returns false, it is still possible to call
% getPrecon, as the default preconditioner is simply the identity operator.
% This check function is mostly useful to tell whether that default
% preconditioner will be in use or not.
%
% See also: canGetDirectionalDerivative canGetGradient canGetHessian

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, July 3, 2013.
% Contributors: 
% Change log: 

    candoit = isfield(problem, 'precon');
    
end
