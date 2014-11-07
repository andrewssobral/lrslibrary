function loc = subdims(subd, t)
%SUBDIMS Compute the locations of subscripts within a subdimension.
%
%  LOC = SUBDIMS(SUBD,T) finds the locations of the subscripts in T
%  that are within the range specified by the cell array SUBD. For
%  example, if SUBD = {1, [1,2], [1,2]}, then the locations of
%  all elements of T that have a first subscript equal to 1, a
%  second subscript equal to 1 or 2, and a third subscript equal to
%  1 or 2 is returned.
%
%MATLAB Tensor Toolbox.
%Copyright 2012, Sandia Corporation.

% This is the MATLAB Tensor Toolbox by T. Kolda, B. Bader, and others.
% http://www.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2012) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in the file LICENSE.txt


% Error check that subd is the right size
if size(subd,2) ~= ndims(t)
    error('Number of subdimensions must equal number of dimensions');
end

% Error check that range is valid
for i = 1 : ndims(t)
    r = subd{i};
    okcolon = ischar(r) && r == ':';
    oknumeric = isreal(r) && ~any(isnan(r(:))) && ~any(isinf(r(:))) ...
        && isequal(r,round(r)) && all(r(:) > 0);        
    if ~(okcolon || oknumeric)
        error('Invalid subdimension.');
    end
end

% Copy out the subscripts from t
subs = t.subs;

if isempty(subs)
    loc = [];
    return;
end

% Compute the indices of the subscripts that are within the
% specified range. We start with all indices in loc and
% pare it down to a final list.

loc = (1:size(subs,1))';
for i = 1:ndims(t)
    if ~(ischar(subd{i}) && subd{i} == ':')

        % Find the subscripts that match in dimension i
        tf = ismember(subs(loc,i), subd{i});

        % Pare down the list of indices
        loc = loc(tf);
    
    end
end

