function C = plus(A,B)
%PLUS Binary addition for sparse tensors. 
%
%   PLUS(A,B) is called for the syntax 'A + B' when A or B is a sparse
%   tensor. A and B must have the same size, unless one is a scalar. A
%   scalar can be added to a sparse tensor of any size.   
%
%   Examples
%   A = sptenrand([4 3 2],5); B = sptenrand([4 3 2],3);
%   A + B %<-- sparse
%   A + 5 %<-- dense
%   A + 0 %<-- dense
%   A + full(A) %<-- dense
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
% The result of a + 5 is dense!
% The result of a + 0 is dense!
% The result of a + full(a) is dense!
% The result of a + b (two sparse matrices) is sparse.

%% Case 1: One argument is a scalar
% Emulating the sparse matrix case here, which creates and returns
% a dense result, even if the scalar is zero.

% Case 1a: Second argument is a scalar or a dense tensor
if isscalar(B) || isa(B,'tensor')
    C = full(A) + B;
    return;
end

% Case 1b: First argument is a scalar
if isscalar(A) 
    C = A + full(B);
    return;
end

%% Case 2: Both are sparse tensors

if ~isa(A,'sptensor') || ~isa(B,'sptensor') || ~isequal(size(A),size(B))
    error('Must be two sparse tensors of the same size');
end

C = sptensor([A.subs; B.subs], [A.vals; B.vals], size(A));
