function [f,g] = example2(x,Data)
%EXAMPLE2   Matrix factorization.
% Compute Frobenius norm residual of an approximate 
% (rank-reduced) two-factor decomposition of a matrix:
%
%     1/2 * ( || A - U*V' ||_F )^2
%
% Derivatives are computed using the matrix form.
%
% Input
%     Data.A:    matrix being approximated, A
%     Data.rank: rank of approximation factors
%     x:         approximation encoded as follows:
%                U = reshape(x(1:m*k),m,k);
%                V = reshape(x(m*k+1:m*k+n*k),n,k);
%
% Output
%     f: function value (residual)
%     g: first derivatives of f w.r.t. x (i.e., U & V)
%
% Matrix dimensions
%     A is m x n
%     U is m x k
%     V is n x k
%
% Examples
%     % Defaults for example2_init: m = 10; n = 8; k = 6;
%     [x0,Data] = example2_init;
%     [f,g] = example2(x0,Data)
%     
%     % Larger problem
%     m = 100; n = 80; k = 10;
%     [x0,Data] = example2_init(m,n,k);
%     out = ncg(@(x) example2(x,Data), x0)
%
%     % Increase amount of computation allowed
%     params = out.Params.Results;
%     params.MaxIters = 1000;
%     params.MaxFuncEvals = 2000;
%     out = ncg(@(x) example2(x,Data), x0, params)
%
%MATLAB Poblano Toolbox.
%Copyright 2009-2012, Sandia Corporation.

% Notes:
%     - Derivatives are computed using matrix form. 
%


%% Data setup
% Data.A should have matrix being modeled
[m,n] = size(Data.A);
k = Data.rank;
U = reshape(x(1:m*k),m,k);
V = reshape(x(m*k+1:m*k+n*k),n,k);

%% Function value (residual)
AmUVt = Data.A-U*V';
f = 0.5*norm(AmUVt,'fro')^2;

%% First derivatives computed in matrix form
g = zeros((m+n)*k,1);
g(1:m*k) = -reshape(AmUVt*V,m*k,1);
g(m*k+1:end) = -reshape(AmUVt'*U,n*k,1);

