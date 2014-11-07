function [sdims,vidx] = tt_dimscheck(dims,N,M)
%TT_DIMSCHECK Used to preprocess dimensions tensor dimensions.
%
%   NEWDIMS = TT_DIMCHECK(DIMS,N) checks that the specified dimensions
%   are valid for a tensor of order N. If DIMS is empty, then
%   NEWDIMS=1:N. If DIMS is negative, then NEWDIMS is everything
%   but the dimensions specified by -DIMS. Finally, NEWDIMS is
%   returned in sorted order.
%
%   [NEWDIMS,IDX] = TT_DIMCHECK(DIMS,N,M) does all of the above but
%   also returns an index for M muliplicands. 
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


% Fix empty case
if isempty(dims)
    dims = 1:N;
end

% Fix "minus" case
if (max(dims) < 0)
    % Check that every member of dims is in 1:N
    tf = ismember(-dims,1:N);
    if min(tf) == 0
        error('Invalid dimensions specified');
    end
    dims = setdiff(1:N, -dims);
end

% Check that every member of dims is in 1:N
tf = ismember(dims,1:N);
if min(tf) == 0
    error('Invalid dimensions specified');
end

% Save the number of dimensions in dims
P = length(dims);

% Reorder dims from smallest to largest (this matters in particular
% for the vector multiplicand case, where the order affects the
% result)
[sdims,sidx] = sort(dims,'ascend');

if (nargout == 2)
    % Can't have more multiplicands them dimensions
    if (M > N)
        error('Cannot have more multiplcands than dimensions');
    end
    
    % Check that the number of mutliplicands must either be
    % full-dimensional (i.e., M==N) or equal to the number of specified
    % dimensions (i.e., M==P).
    if (M ~= N) && (M ~= P)
        error('Invalid number of multiplicands');
    end
    
    % Check sizes to determine how to index multiplicands
    if (P == M)
        % Case 1: Number of items in dims and number of multiplicands
        % are equal; therefore, index in order of how sdims was sorted.
        vidx = sidx;   
    else
        % Case 2: Number of multiplicands is equal to the number of
        % dimensions in the tensor; therefore, index multiplicands by
        % dimensions specified in dims argument.
        vidx = sdims;  % index multiplicands by (sorted) dimension
    end
end
