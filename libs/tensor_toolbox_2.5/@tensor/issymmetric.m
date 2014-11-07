function [tf,all_diffs,all_perms] = issymmetric(X,grps)
%ISSYMMETRIC Verify that a tensor X is symmetric in specified modes.
%
%   TF = ISSYMMETRIC(X) returns true if X is exactly symmetric for every
%   permutation.
%
%   [TF,DIFFS,PERMS] = ISSYMMETRIC(X) also returns that maximum difference
%   in DIFFS for each permutation in PERMS (one permutation per row).
%
%   [...] = ISSYMMETRIC(X,IDX) checks symmetry with respect to the modes
%   specified in IDX, which can be an array of indices or a cell array of
%   arrays of symmetric indices.
%
%   See also SYMMETRIZE.
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

%T. Kolda, April 2011.

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

% Check tensor dimensions for compatibility with symmetrization
for i = 1:length(grps)
    dims = grps{i};
    for j = dims(2:end)
        if sz(j) ~= sz(dims(1))
            tf = false; 
            return;
        end
    end
end

% Check actual symmetry
cnt = sum(cellfun(@(x) factorial(length(x)), grps));
all_diffs = zeros(cnt,1);
all_perms = zeros(cnt,n);
idx = 1;
for i = 1:length(grps)
    
    % Compute the permutations for this group of symmetries
    p = perms(grps{i});
    
    for j = 1:size(p,1)
        
        % Create the permutation to check
        q = 1:n;
        q(grps{i}) = p(j,:);
        
        % Save the permutation        
        all_perms(idx,:) = q;       
        
        % Do the permutation and see if it's a match. If it's not a match,
        % record the difference.
        Y = permute(X,q);
        if isequal(X.data,Y.data)
            all_diffs(idx) = 0;
        else
            all_diffs(idx) = max(abs(X.data(:)-Y.data(:)));
        end
    
        % Increment the index
        idx = idx + 1;
        
    end
    
end

tf = all(all_diffs == 0);
