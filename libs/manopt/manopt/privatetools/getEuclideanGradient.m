function [egrad, storedb] = getEuclideanGradient(problem, x, storedb)
% Computes the Euclidean gradient of the cost function at x.
%
% function [egrad, storedb] = getEuclideanGradient(problem, x, storedb)
%
% Returns the Euclidean gradient at x of the cost function described in the
% problem structure. The cache database storedb is passed along, possibly
% modified and returned in the process.
%
% Because computing the Hessian based on the Euclidean Hessian will require
% the Euclidean gradient every time, to avoid overly redundant
% computations, if the egrad function does not use the store caching
% capabilites, we implement an automatic caching functionality. This means
% that even if the user does not use the store parameter, the hashing
% function will be used, and this can translate in a performance hit for
% small problems. For problems with expensive cost functions, this should
% be a bonus though. Writing egrad to accept the optional store parameter
% (as input and output) will disable automatic caching, but activate user
% controlled caching, which means hashing will be computed in all cases.
%
% If you absolutely do not want hashing to be used (and hence do not want
% caching to be used), you can define grad instead of egrad, without store
% support, and call problem.M.egrad2rgrad manually.
%
% See also: getGradient canGetGradient canGetEuclideanGradient

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, July 9, 2013.
% Contributors: 
% Change log: 

    
    if isfield(problem, 'egrad')
    %% Compute the Euclidean gradient using egrad.

		is_octave = exist('OCTAVE_VERSION', 'builtin');
		if ~is_octave
			narg = nargin(problem.egrad);
		else
			narg = 2;
		end
	
        % Check whether the egrad function wants to deal with the store
        % structure or not.
        switch narg
            case 1
                % If it does not want to deal with the store structure,
                % then we do some caching of our own. There is a small
                % performance hit for this is some cases, but we expect
                % that this is most often the preferred choice.
                store = getStore(problem, x, storedb);
                if ~isfield(store, 'egrad__')
                    store.egrad__ = problem.egrad(x);
                    storedb = setStore(problem, x, storedb, store);
                end
                egrad = store.egrad__;
            case 2
                % Obtain, pass along, and save the store structure
                % associated to this point. If the user deals with the
                % store structure, then we don't do any automatic caching:
                % the user is in control.
                store = getStore(problem, x, storedb);
                [egrad, store] = problem.egrad(x, store);
                storedb = setStore(problem, x, storedb, store);
            otherwise
                up = MException('manopt:getEuclideanGradient:badegrad', ...
                    'egrad should accept 1 or 2 inputs.');
                throw(up);
        end

    else
    %% Abandon computing the Euclidean gradient
    
        up = MException('manopt:getEuclideanGradient:fail', ...
            ['The problem description is not explicit enough to ' ...
             'compute the Euclidean gradient of the cost.']);
        throw(up);
        
    end
    
end
