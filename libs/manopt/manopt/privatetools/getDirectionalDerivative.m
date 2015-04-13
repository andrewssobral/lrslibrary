function [diff, storedb] = getDirectionalDerivative(problem, x, d, storedb)
% Computes the directional derivative of the cost function at x along d.
%
% function [diff, storedb] = getDirectionalDerivative(problem, x, d, storedb)
%
% Returns the derivative at x along d of the cost function described in the
% problem structure. The cache database storedb is passed along, possibly
% modified and returned in the process.
%
% See also: getGradient canGetDirectionalDerivative

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 

    
    if isfield(problem, 'diff')
    %% Compute the directional derivative using diff.
        
		is_octave = exist('OCTAVE_VERSION', 'builtin');
		if ~is_octave
			narg = nargin(problem.diff);
		else
			narg = 3;
		end
		
        % Check whether the diff function wants to deal with the store
        % structure or not.
        switch narg
            case 2
                diff = problem.diff(x, d);
            case 3
                % Obtain, pass along, and save the store structure
                % associated to this point.
                store = getStore(problem, x, storedb);
                [diff store] = problem.diff(x, d, store);
                storedb = setStore(problem, x, storedb, store);
            otherwise
                up = MException('manopt:getDirectionalDerivative:baddiff', ...
                    'diff should accept 2 or 3 inputs.');
                throw(up);
        end
    
    elseif canGetGradient(problem)
    %% Compute the directional derivative using the gradient.
        
        % Compute the gradient at x, then compute its inner product with d.
        [grad, storedb] = getGradient(problem, x, storedb);
        diff = problem.M.inner(x, grad, d);
        
    else
    %% Abandon computing the directional derivative.
    
        up = MException('manopt:getDirectionalDerivative:fail', ...
            ['The problem description is not explicit enough to ' ...
             'compute the directional derivatives of f.']);
        throw(up);
        
    end
    
end
