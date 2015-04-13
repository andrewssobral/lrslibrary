function stats = applyStatsfun(problem, x, storedb, options, stats)
% Apply the statsfun function to a stats structure (for solvers).
%
% function stats = applyStatsfun(problem, x, storedb, options, stats)
%
% Applies the options.statsfun user supplied function to the stats
% structure, if it was provided, with the appropriate inputs, and returns
% the (possibly) modified stats structure.
%
% See also: 

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, April 3, 2013.
% Contributors: 
% Change log: 

	if isfield(options, 'statsfun')

		is_octave = exist('OCTAVE_VERSION', 'builtin');
		if ~is_octave
			narg = nargin(options.statsfun);
		else
			narg = 4;
		end
		
        switch narg
            case 3
                stats = options.statsfun(problem, x, stats);
            case 4
                store = getStore(problem, x, storedb);
                stats = options.statsfun(problem, x, stats, store);
            otherwise
                warning('manopt:statsfun', ...
                        'statsfun unused: wrong number of inputs');
        end
	end

end
