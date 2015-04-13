function candoit = canGetEuclideanGradient(problem)
% Checks whether the Euclidean gradient can be computed for a problem.
%
% function candoit = canGetEuclideanGradient(problem)
%
% Returns true if the Euclidean gradient can be computed given the problem
% description, false otherwise.
%
% See also: canGetGradient getEuclideanGradient

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 


    candoit = isfield(problem, 'egrad');
    
end
