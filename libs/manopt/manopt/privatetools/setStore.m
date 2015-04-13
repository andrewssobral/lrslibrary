function storedb = setStore(problem, x, storedb, store)
% Updates the store struct. pertaining to a point in the storedb database.
%
% function storedb = setStore(problem, x, storedb, store)
%
% Updates the storedb database of structures such that the structure
% corresponding to the point x will be replaced by store. If there was no
% record for the point x, it is created and set to store. The updated
% storedb database is returned. The lastset__ field of the store structure
% keeps track of which stores were updated latest.

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 
%   Dec. 6, 2013, NB:
%       Now using a persistent uint32 counter instead of cputime to track
%       the most recently modified stores.

    % This persistent counter is used to keep track of the order in which
    % store structures are updated. This is used by purgeStoredb to erase
    % the least recently useful store structures first when garbage
    % collecting.
    persistent counter;
    if isempty(counter)
        counter = uint32(0);
    end

    assert(nargout == 1, ...
           'The output of setStore should replace your storedb.');
   
    % Construct the fieldname (key) associated to the current point x.
    key = problem.M.hash(x);
    
    % Set the value associated to that key to store.
    storedb.(key) = store;
    
    % Add / update a last-set flag.
    storedb.(key).lastset__ = counter;
    counter = counter + 1;

end
