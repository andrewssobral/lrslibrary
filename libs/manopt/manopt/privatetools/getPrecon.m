function [Pd, storedb] = getPrecon(problem, x, d, storedb)
% Applies the preconditioner for the Hessian of the cost at x along d.
%
% function [Pd, storedb] = getPrecon(problem, x, storedb)
%
% Returns as Pd the result of applying the Hessian preconditioner to the
% tangent vector d at point x. If no preconditioner is specified, Pd = d
% (identity).
%
% See also: getHessian

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 

    
    if isfield(problem, 'precon')
    %% Compute the preconditioning using precon.

		is_octave = exist('OCTAVE_VERSION', 'builtin');
		if ~is_octave
			narg = nargin(problem.precon);
		else
			narg = 3;
		end
	
        % Check whether the precon function wants to deal with the store
        % structure or not.
        switch narg
            case 2
                Pd = problem.precon(x, d);
            case 3
                % Obtain, pass along, and save the store structure
                % associated to this point.
                store = getStore(problem, x, storedb);
                [Pd store] = problem.precon(x, d, store);
                storedb = setStore(problem, x, storedb, store);
            otherwise
                up = MException('manopt:getPrecon:badprecon', ...
                    'precon should accept 2 or 3 inputs.');
                throw(up);
        end      

    else
    %% No preconditioner provided, so just use the identity.
    
        Pd = d;
        
    end
    
end
