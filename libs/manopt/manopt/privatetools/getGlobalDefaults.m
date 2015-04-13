function opts = getGlobalDefaults()
% Returns a structure with default option values for Manopt.
%
% function opts = getGlobalDefaults()
%
% Returns a structure opts containing the global default options such as
% verbosity level etc. Typically, global defaults are overwritten by solver
% defaults, which are in turn overwritten by user-specified options.
% See the online Manopt documentation for details on options.
%
% See also: mergeOptions

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 


    % Verbosity level: 0 is no output at all. The higher the verbosity, the
    % more info is printed / displayed during solver execution.
    opts.verbosity = 3;
    
    % If debug is set to true, additional computations may be performed and
    % debugging information is outputed during solver execution.
    opts.debug = false;
    
    % Maximum number of store structures to store. If set to 0, caching
    % capabilities are not disabled, but the cache will be emptied at each
    % iteration of iterative solvers (more specifically: every time the
    % solver calls the purgeStoredb tool).
    opts.storedepth = 20;
    
    % Maximum amount of time a solver may execute, in seconds.
    opts.maxtime = inf;

end
