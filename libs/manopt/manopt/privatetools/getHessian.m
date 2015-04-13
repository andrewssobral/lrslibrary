function [hess, storedb] = getHessian(problem, x, d, storedb)
% Computes the Hessian of the cost function at x along d.
%
% function [hess, storedb] = getHessian(problem, x, d, storedb)
%
% Returns the Hessian at x along d of the cost function described in the
% problem structure. The cache database storedb is passed along, possibly
% modified and returned in the process.
%
% If an exact Hessian is not provided, an approximate Hessian is returned
% if possible, without warning. If not possible, an exception will be
% thrown. To check whether an exact Hessian is available or not (typically
% to issue a warning if not), use canGetHessian.
%
% See also: getPrecon getApproxHessian canGetHessian

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 
    
    if isfield(problem, 'hess')
    %% Compute the Hessian using hess.

		is_octave = exist('OCTAVE_VERSION', 'builtin');
		if ~is_octave
			narg = nargin(problem.hess);
		else
			narg = 3;
		end
	
        % Check whether the hess function wants to deal with the store
        % structure or not.
        switch narg
            case 2
                hess = problem.hess(x, d);
            case 3
                % Obtain, pass along, and save the store structure
                % associated to this point.
                store = getStore(problem, x, storedb);
                [hess, store] = problem.hess(x, d, store);
                storedb = setStore(problem, x, storedb, store);
            otherwise
                up = MException('manopt:getHessian:badhess', ...
                    'hess should accept 2 or 3 inputs.');
                throw(up);
        end
    
    elseif isfield(problem, 'ehess') && canGetEuclideanGradient(problem)
    %% Compute the Hessian using ehess.
    
        % We will need the Euclidean gradient for the conversion from the
        % Euclidean Hessian to the Riemannian Hessian.
        [egrad, storedb] = getEuclideanGradient(problem, x, storedb);
        
		is_octave = exist('OCTAVE_VERSION', 'builtin');
		if ~is_octave
			narg = nargin(problem.ehess);
		else
			narg = 3;
		end
		
        % Check whether the ehess function wants to deal with the store
        % structure or not.
        switch narg
            case 2
                ehess = problem.ehess(x, d);
            case 3
                % Obtain, pass along, and save the store structure
                % associated to this point.
                store = getStore(problem, x, storedb);
                [ehess, store] = problem.ehess(x, d, store);
                storedb = setStore(problem, x, storedb, store);
            otherwise
                up = MException('manopt:getHessian:badehess', ...
                    'ehess should accept 2 or 3 inputs.');
                throw(up);
        end
        
        % Convert to the Riemannian Hessian
        hess = problem.M.ehess2rhess(x, egrad, ehess, d);
        
    else
    %% Attempt the computation of an approximation of the Hessian.
        
        [hess, storedb] = getApproxHessian(problem, x, d, storedb);
        
    end
    
end
