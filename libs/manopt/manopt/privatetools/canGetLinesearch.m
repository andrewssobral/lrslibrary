function candoit = canGetLinesearch(problem)
% Checks whether the problem structure can give a line-search a hint.
%
% function candoit = canGetLinesearch(problem)
%
% Returns true if the the problem description includes a mechanism to give
% line-search algorithms a hint as to "how far to look", false otherwise.
%
% See also: getLinesearch

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, July 17, 2014.
% Contributors: 
% Change log: 


    candoit = isfield(problem, 'linesearch');
    
end
