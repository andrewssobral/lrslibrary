function [approxhess, storedb] = getApproxHessian(problem, x, d, storedb)
% Computes an approximation of the Hessian of the cost fun. at x along d.
%
% function [approxhess, storedb] = getApproxHessian(problem, x, d, storedb)
%
% Returns an approximation of the Hessian at x along d of the cost function
% described in the problem structure. The cache database storedb is passed
% along, possibly modified and returned in the process.
%
% If no approximate Hessian was furnished, this call is redirected to
% getHessianFD.
%
% See also: getHessianFD

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 


    if isfield(problem, 'approxhess')
    %% Compute the approximate Hessian using approxhess.
        
		is_octave = exist('OCTAVE_VERSION', 'builtin');
		if ~is_octave
			narg = nargin(problem.approxhess);
		else
			narg = 3;
		end
		
        % Check whether the approximate Hessian function wants to deal with
        % the store structure or not.
        switch narg
            case 2
                approxhess = problem.approxhess(x, d);
            case 3
                % Obtain, pass along, and save the store structure
                % associated to this point.
                store = getStore(problem, x, storedb);
                [approxhess store] = problem.approxhess(x, d, store);
                storedb = setStore(problem, x, storedb, store);
            otherwise
                up = MException('manopt:getApproxHessian:badapproxhess', ...
                    'approxhess should accept 2 or 3 inputs.');
                throw(up);
        end
        
    else
    %% Try to fall back to a standard FD approximation.
    
        [approxhess, storedb] = getHessianFD(problem, x, d, storedb);
        
    end
    
end
