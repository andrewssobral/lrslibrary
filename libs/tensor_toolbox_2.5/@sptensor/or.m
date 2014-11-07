function C = or(A,B)
%OR Logical OR (|) for sptensors.
%
%   See also SPTENSOR.
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


%% Observations for sparse matrix case.
% The result of a | 5 is dense!
% The result of a | 0 is dense!
% The result of a | full(a) is dense!
% The result of a | a is sparse.

%% Case 1: One argument is a scalar
if isscalar(B) || isa(B,'tensor')
    C = full(A) | B;
    return;
end
if isscalar(A)
    C = A | full(B);
    return;
end

%% Case 2: Both A and B are sparse tensors
if ~isequal(size(A),size(B))
    error('Must be tensors of the same size');
end

if isa(A,'sptensor') && isa(B,'sptensor')
    C = sptensor([A.subs; B.subs], 1, size(A), @(x) length(x) >= 1);
    return;
end

%% Otherwise
error('The arguments must be two sptensors or an sptensor and a scalar.');
