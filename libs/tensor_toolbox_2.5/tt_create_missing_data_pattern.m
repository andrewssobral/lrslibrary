function W = create_missing_data_pattern(sz,M,isSparse)
%TEST_CREATE_RME Creates a randomly missing element (RME) indicator tensor.
%
%   W = TEST_CREATE_RME(SZ,M) creates an indicator (binary) tensor W of the
%   specified size with 0's indicating missing data and 1's indicating
%   valid data. The percentage of zeros is given by M. Will only return a
%   tensor that has at least one entry per N-1 dimensional slice. 
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

%   Code by Evrim Acar and Tammy Kolda, 2009.

%% Set up isSparse variable
if ~exist('isSparse','var')
    isSparse = false;
end

%% Initialize
% Number of dimensions
N = length(sz);

% Total number of entries in tensor of given sz
P = prod(sz);

% Total number of entries that should be set to one
Q = ceil((1-M)*P);

%% Create the tensor
% Keep iterating until the tensor is created or we give up.
for iter = 1:20
    % Create the indicator tensor W
    if isSparse
        % start with 50% more than Q random subs
        % TODO: work out the expected value of a*Q to guarantee Q unique entries
        subs = unique(ceil(rand(ceil(1.5*Q),size(sz,2))*diag(sz)),'rows');
        % check if there are too many unique subs
        if size(subs,1) > Q
            % unique orders the subs and would bias toward first subs
            % with lower values, so we sample to cut back
            idx = randperm(size(subs,1));
            subs = subs(idx(1:Q),:);
        elseif size(subs,1) < Q
            warning('Only generated %d of %d desired subscripts', size(subs,1), Q);
        end
        W = sptensor(subs,1,sz);
    else
        % Compute the linear indices of the missing entries. Note that
        % the indices must be a column array for the linear indexing
        % into W to  work.
        idx = randperm(P);
        idx = idx(1:Q)';
        W = tenzeros(sz);
        W(idx) = 1;
    end
    
    % Check if W has any empty slices
    isokay = zeros(N,1);
    for n = 1:N
        isokay(n) = all(double(collapse(W,-n)));
    end 

    % Quit if we're okay
    if all(isokay)
         break;
    end
        
end

if ~all(isokay)
    error('After %d iterations, cannot produce a tensor with %f%% missing data without an empty slice.', iter, M*100);
end

    