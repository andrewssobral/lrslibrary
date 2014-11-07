function x = fac_to_vec(A)
%FAC_TO_VEC Converts a set of factor matrices to a vector.
%
%   X = FAC_TO_VEC(A) converts a cell array of factor matrices A to a
%   vector by vectorizing each matrix and stacking them.
%
%   See also CP_VEC_TO_FAC, CP_FUN, CP_OPT.
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


%% Set-up
N = length(A);

%% Get sizes
sz = zeros(N,1);
for n = 1:N
    sz(n) = size(A{n},1);
end
R = size(A{1},2);
P = sum(sz)*R;

%% Create x
x = zeros(P,1);
for n = 1:N
    idx1 = sum(sz(1:n-1))*R + 1;
    idx2 = sum(sz(1:n))*R;
    x(idx1:idx2) = reshape(A{n},sz(n)*R,1);
end