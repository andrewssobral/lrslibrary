function Y = symmetrize(X,grps)
%SYMMETRIZE Symmetrize a tensor X in specified modes.
%
%   Y = symmetrize(X) will symmetrize a tensor X with respect to all
%   modes so that Y is symmetric with respect to any permutation of
%   indices.
%
%   Y = symmetrize(X,MODES) will symmetrix a tensor X with respect to the
%   modes specified by the vector MODES of mode indices. The second
%   argument may alternatively be a cell array of vectors of modes to,
%   e.g., specify that it should be symmetric with respect to mode [1 3] as
%   well as [2 4].
%
%   See also ISSYMMETRIC.
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

%T. Kolda, April 2011
n = ndims(X);
sz = size(X);

% Check that grps exists; if not, create it.
if ~exist('grps','var')
    grps = {1:n};
end

% Check that grps is a cell array.
if ~iscell(grps)
    grps = {grps};
end

if ~isnumeric(grps{1})
    error('MODES must be numeric');
end

% Check tensor dimensions for compatibility with symmetrization
ngrps = length(grps);
for i = 1:ngrps
    dims = grps{i};
    for j = dims(2:end)
        if sz(j) ~= sz(dims(1))
            error('Dimension mismatch for symmetrization');
        end
    end
end

% Check for no overlap in the sets
for i = 1:ngrps
    for j = i+1:ngrps
        if ~isempty(intersect(grps{i},grps{j}))
            error('Cannot haver overlapping symmetries');
        end
    end
end

% Create the combinations for each symmetrized subset
combos = cell(ngrps,1);
for i = 1:ngrps
    combos{i} = perms(grps{i});   
end

% Create all the permuations to be averaged
total_perms = prod(cellfun(@length,combos));
sym_perms = repmat(1:n, total_perms, 1);
for i = 1:ngrps
    ntimes = prod(cellfun(@length,combos(1:i-1))); 
    ncopies = prod(cellfun(@length,combos(i+1:end))); 
    nelems = length(combos{i});
    
    idx = 1;
    for j = 1:ntimes
        for k = 1:nelems
            for l = 1:ncopies
                sym_perms(idx,grps{i}) = combos{i}(k,:);
                idx = idx + 1;
            end
        end
    end
end

% Create an average tensor
Y = tenzeros(size(X));
for i = 1:total_perms
    Y = Y + permute(X,sym_perms(i,:));    
end
Y = Y / total_perms;

% It's not *exactly* symmetric due to oddities in differently ordered
% summations and so on, so let's fix that.
% Idea borrowed from Gergana Bounova:
% http://www.mit.edu/~gerganaa/downloads/matlab/symmetrize.m
for i = 1:total_perms
    Z = permute(Y,sym_perms(i,:));
    Y.data(:) = max(Y.data(:),Z.data(:));    
end



    



