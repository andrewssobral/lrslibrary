function [D, alpha, E, A] = inexact_alm_rosl_subsampling(X, K, L, lambda, tol, maxIter) 

% This matlab code implements the inexact augmented Lagrange multiplier 
% method for Robust Orthogonal Subspace Learning with random sampling.
%
% X - m x n matrix of observations/data (required input)
%
% lambda - weight on sparse error term in the cost function
%
% tol - tolerance for stopping criterion.
%
% maxIter - maximum number of iterations
%  
%  
% min \sum_{1<i<T}\|\alpha_i\|_2 + \lambda\|E\|_1, s.t. \|D_i\|_2= 1 ;
% 
% The laplacian function: 
% L(A,E,Y,u) =  \sum_1^T |alpha_i|_2 + \lambda*|E|_1 + \mu/2* |X- \sum_1^T D_i*alpha_i-E+Y/mu|_F^2;
%
% Xianbiao Shu (xshu2@illinois.edu)
% Copyright: Mitsubishi Electric Research Lab
% Reference: X. Shu, F. Porikli, N. Ahujia "Robust Orthonormal Subspace Learning: Efficient Recovery of Corrupted Low-rank Matrices". CVPR 2014; 

[m, n] = size(X);


% Random subsampling (using randsample or randperm)
 
row_sampled = randsample(m,L);
row_notsampled = [1:m]';
row_notsampled(row_sampled) =[];
row_ind = [row_sampled; row_notsampled];
column_sampled = randsample(n,L);
column_notsampled = [1:n]';
column_notsampled(column_sampled) =[];
column_ind = [column_sampled; column_notsampled];

X_perm = X(row_ind, :); 
X_perm = X_perm(:, column_ind);



% % solve the simplified RSL
% [Dict_hat1, alpha] = inexact_alm_rsl(X_perm(1:L,:), K, lambda, tol, maxIter);
% 
% % robust linear regression
% [Dict_hat2] = inexact_alm_rlr(X_perm(L+1:m,1:L), K, alpha(:, 1:L), tol, maxIter) ;


% solve the simplified RSL
[D, alpha_1] = inexact_alm_rosl(X_perm(:,1:L), K, lambda, tol, maxIter);

% robust linear regression
[alpha] = inexact_alm_rlr(X_perm(1:L, :), K, D(1:L, :), tol, maxIter) ;

% compute the low-rank matrix 
%  D = [Dict_hat1; Dict_hat2]; 
 A = D*alpha;
 
 %permute back
 A(:, column_ind) = A; 
 A(row_ind, :) = A;
 
 E = X -A; 
 
 


end


