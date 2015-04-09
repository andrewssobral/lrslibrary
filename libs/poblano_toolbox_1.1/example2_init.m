function [x0,Data] = example_init(m,n,k,A)
%EXAMPLE2_INIT   Matrix factorization example: initialization.
%
% Initialize data for use in computing an approximate 
% (rank-reduced) two-factor decomposition of a matrix:
%
%     A \approx U*V'
%
% Input
%     m: number of rows of matrix to be approximated (10)
%     n: number of columns of matrix to be approximated (8)
%     k: rank of approximation factors (6)
%     A: matrix to be approximated (optional)
%
% Output
%     x0:   vector of variables in two-factor approximation 
%     Data: A: matrix to be approximated
%           rank: rank of approximation factors (== k)
%
% Matrix dimensions
%     Data.A is m x n
%     x0     is (m+n)*k x 1
%
% Example
%     m = 3; n = 4; k = 2;
%     [x0,Data] = example2_init(m,n,k)
%
%MATLAB Poblano Toolbox.
%Copyright 2009-2012, Sandia Corporation.

% Notes:
%      - Ordering of vector of variables must be consistent across uses
%     

%% Check parameters
if (nargin < 3)
    k = 6;
end
if (nargin < 2)
    n = 8;
end
if (nargin < 1)
    m = 10;
end

%% Check to see if dimensions are valid
if (k > min(m,n))
    error('The value of k is too large for the values of m and n!');
end

%% Data matrix and rank of approximation
Data.rank = k;
if (nargin < 4)
    Data.A = randn(m,n);
elseif ( (m ~= size(A,1)) && (n ~= size(A,2)) )
    error('A is not an m x n matrix');
else
    Data.A = A;
end

%% Initial two-factor approximation
U0 = randn(m,k);
V0 = randn(n,k);
U0 = reshape(U0,m*k,1);
V0 = reshape(V0,n*k,1);
x0 = [U0;V0];
