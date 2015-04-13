function [cost, storedb] = getCost(problem, x, storedb)
% Computes the cost function at x.
%
% function [cost, storedb] = getCost(problem, x, storedb)
%
% Returns the value at x of the cost function described in the problem
% structure. The cache database storedb is passed along, possibly modified
% and returned in the process.
%
% See also: canGetCost

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 


    if isfield(problem, 'cost')
    %% Compute the cost function using cost.

		is_octave = exist('OCTAVE_VERSION', 'builtin');
		if ~is_octave
			narg = nargin(problem.cost);
		else
			narg = 2;
		end
	
        % Check whether the cost function wants to deal with the store
        % structure or not.
        switch narg
            case 1
                cost = problem.cost(x);
            case 2
                % Obtain, pass along, and save the store structure
                % associated to this point.
                store = getStore(problem, x, storedb);
                [cost, store] = problem.cost(x, store);
                storedb = setStore(problem, x, storedb, store);
            otherwise
                up = MException('manopt:getCost:badcost', ...
                    'cost should accept 1 or 2 inputs.');
                throw(up);
        end
        
    elseif isfield(problem, 'costgrad')
    %% Compute the cost function using costgrad.

		is_octave = exist('OCTAVE_VERSION', 'builtin');
		if ~is_octave
			narg = nargin(problem.costgrad);
		else
			narg = 2;
		end
	
        % Check whether the costgrad function wants to deal with the store
        % structure or not.
        switch narg
            case 1
                cost = problem.costgrad(x);
            case 2
                % Obtain, pass along, and save the store structure
                % associated to this point.
                store = getStore(problem, x, storedb);
                [cost, grad, store] = problem.costgrad(x, store); %#ok
                storedb = setStore(problem, x, storedb, store);
            otherwise
                up = MException('manopt:getCost:badcostgrad', ...
                    'costgrad should accept 1 or 2 inputs.');
                throw(up);
        end

    else
    %% Abandon computing the cost function.

        up = MException('manopt:getCost:fail', ...
            ['The problem description is not explicit enough to ' ...
             'compute the cost.']);
        throw(up);
        
    end
    
end
