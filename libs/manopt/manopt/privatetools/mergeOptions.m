function opts = mergeOptions(opts1, opts2)
% Merges two options structures with one having precedence over the other.
%
% function opts = mergeOptions(opts1, opts2)
%
% input: opts1 and opts2 are two structures.
% output: opts is a structure containing all fields of opts1 and opts2.
% Whenever a field is present in both opts1 and opts2, it is the value in
% opts2 that is kept.
%
% The typical usage is to have opts1 contain default options and opts2
% contain user-specified options that overwrite the defaults.
%
% See also: getGlobalDefaults

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 


    if isempty(opts1)
        opts1 = struct();
    end
    if isempty(opts2)
        opts2 = struct();
    end

    opts = opts1;
    fields = fieldnames(opts2);
    for i = 1 : length(fields)
        opts.(fields{i}) = opts2.(fields{i});
    end
    
end
