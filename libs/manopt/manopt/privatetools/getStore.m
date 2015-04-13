function store = getStore(problem, x, storedb)
% Extracts a store struct. pertaining to a point from the storedb database.
%
% function store = getStore(problem, x, storedb)
%
% Queries the storedb database of structures (itself a structure) and
% returns the store structure corresponding to the point x. If there is no
% record for the point x, returns an empty structure.
%
% See also: setStore purgeStoredb

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 
   
    % Construct the fieldname (key) associated to the queried point x.
    key = problem.M.hash(x);
    
    % If there is a value stored for this key, return it.
    % Otherwise, return an empty structure.
    if isfield(storedb, key)
        store = storedb.(key);
    else
        store = struct();
    end

end
