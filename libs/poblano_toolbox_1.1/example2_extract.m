function [U,V] = example2_extract(m,n,k,x)
%EXAMPLE2_EXTRACT   Matrix factorization example: extraction.
%
% Helper function to extract factors of an approximate 
% (rank-reduced) two-factor decomposition of a matrix:
%
%     A \approx U*V'
%
% Input
%     m: number of rows of U to be approximated
%     n: number of rows of V to be approximated
%     k: rank of approximation factors
%     x: vector of variables in two-factor approximation
%
% Output
%     U: first factor
%     V: second factor
%
% Matrix dimensions
%     A is m x n
%     U is m x k
%     V is n x k
%
% Example
%     m = 3; n = 2; k = 2;
%     [x,Data] = example2_init(m,n,k);
%     [U,V] = example_extract(m,n,k,x)
%     norm(Data.A-U*V')
%
%MATLAB Poblano Toolbox.
%Copyright 2009-2012, Sandia Corporation.

%% Perform check on sizes of inputs
[x_m, x_n] = size(x);
if ( (x_n ~= 1) || (x_m ~= (m+n)*k) )
    error('Dimensions do not agree!');
end

%% Extract approximation factors from x
U = reshape(x(1:m*k),m,k);
V = reshape(x(m*k+1:m*k+n*k),n,k);
