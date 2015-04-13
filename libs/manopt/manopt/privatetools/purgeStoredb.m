function storedb = purgeStoredb(storedb, storedepth)
% Makes sure the storedb database does not exceed some maximum size.
%
% function storedb = purgeStoredb(storedb, storedepth)
%
% Trim the store database storedb such that it contains at most storedepth
% elements (store structures). The 'lastset__' field of the store
% structures is used to delete the oldest elements first.

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 
%   Dec. 6, 2013, NB:
%       Now using a persistent uint32 counter instead of cputime to track
%       the most recently modified stores.


    if storedepth <= 0
        storedb = struct();
        return;
    end

    % Get list of field names (keys).
    keys = fieldnames(storedb);
    nkeys = length(keys);
    
    % If we need to remove some of the elements in the database.
    if nkeys > storedepth
        
        % Get the last-set counter of each element: a higher number means
        % it was modified more recently.
        lastset = zeros(nkeys, 1, 'uint32');
        for i = 1 : nkeys
            store = storedb.(keys{i});
            lastset(i) = store.lastset__;
        end
        
        % Sort the counters and determine the threshold above which the
        % field needs to be removed.
        sortlastset = sort(lastset, 1, 'descend');
        minlastset = sortlastset(storedepth);
        
        % Remove all fields that are too old.
        storedb = rmfield(storedb, keys(lastset < minlastset));
        
    end
    
end
