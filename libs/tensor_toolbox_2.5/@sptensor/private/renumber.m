function [newsubs, newsz] = renumber(subs, sz, range)
%RENUMBER indices for sptensor subsref
%
%  [NEWSUBS,NEWSZ] = RENUMBER(SUBS,SZ,RANGE) takes a set of
%  original subscripts SUBS with entries from a tensor of size
%  SZ. All the entries in SUBS are assumed to be within the
%  specified RANGE. These subscripts are then renumbered so that,
%  in dimension i, the numbers range from 1:numel(RANGE(i)).
%
%  See also SPTENSOR/SUBSREF
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


newsz = sz;
newsubs = subs;
for i = 1 : size(sz,2)
    if ~(ischar(range{i}) && range{i} == ':')
	if (isempty(subs))
	    newsz(i) = numel(range{i});
	else
	    [newsubs(:,i), newsz(i)] = ...
		renumberdim(subs(:,i), sz(i), range{i});
	end
    end
end
	
%------------------------------------------------------
function [newidx, newsz] = renumberdim(idx, sz, range)
%RENUMBERDIM helper function for RENUMBER
%  See also SPTENSOR/PRIVATE/RENUMBER

% Determine the size of the new range
newsz = numel(range);

% Create a map from the old range to the new range
map = zeros(1, sz);
for i = 1 : newsz
    map(range(i)) = i;
end

% Do the mapping
newidx = map(idx);
