function [D, alpha, E, A] = inexact_alm_rosl(X, K, lambda, tol, maxIter) 

% This matlab code implements the inexact augmented Lagrange multiplier 
% method for Robust Orthogonal Subspace Learning.
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

%addpath PROPACK;


[m, n] = size(X);

% initialize
alpha = rand(K,n);  % Low-rank coeffs
D = zeros(m,K);     % Low-rank dictionaries
E = zeros(size(X)); % sparse innovation
A = X;

Z = 0; % Z = Y/mu

norm_inf = norm(X(:), inf); 
d_norm = norm(X, 'fro');

mu =  10*lambda/norm_inf;    %10*lambda/norm_inf; % this parameter can be tuned
rho = 1.5;       %1.3    % this parameter can be tuned
mu_bar = mu * 1e7;

iter = 0;
Iteration = 0;
converged = false;

while ~converged    

    iter = iter + 1;
    
    E_temp = X - A + Z; 

    
    
    %% 3D Total variation threshold using Periodic conditions 
%       E = TV3D_threshold_neumann(E_temp, lambda/mu, h, w);
      
    %% 3D Total variation threshold using Periodic conditions 
%     E = TV3D_threshold_periodic(E_temp, lambda/mu, h, w);
%     E = TV3D_threshold_periodic_zeroPadding(E_temp, lambda/mu, h, w);
      
    %% Total variation threshold using Periodic conditions (good)
%     E = TV_threshold_periodic(E_temp, lambda/mu, h, w);
%     E = TV_threshold_periodic_zeroPadding(E_temp, lambda/mu, h, w);
    
    %% Total variation threshold using Dirichlet conditions (bad)
%     E = TV_threshold_dirichlet(E_temp, lambda/mu, h, w);
%     E = TV_threshold_dirichlet_zeroPadding(E_temp, lambda/mu, h, w);
   
    %% Total variation threshold using Neumann conditions (good)
%     E = TV_threshold_neumann(E_temp, lambda/mu, h, w);   
%     E = TV_threshold_neumann_zeropadding(E_temp, lambda/mu, h, w);


    %% Intensity threshold
    E = max(abs(E_temp) - lambda/mu, 0).*sign(E_temp);   




    
    %% Sparse Representation
    A_ini = X -E + Z;
    
    
%%    min \sum_{1<i<T}\|\alpha_i\|_2+ \lambda*\|E\|_1, s.t. \|D_i\|_2= 1 ;           the best model
    [A, D, alpha] = LowRankDictionaryShrinkage_m(A_ini, D, alpha, 1/mu); 
%     [A, D, alpha] = LowRankDictionaryShrinkage(A_ini, D, alpha, 1/mu); 
    
%%    min (\|D\|_f^2+\sum_{1<i<T}\|\alpha_i\|_2)+ \lambda*\|E\|_1
%     [A, D, alpha] = LowRankDictionaryShrinkage1_m(A_ini, D, alpha, 1/mu); 

%%    min \sum_{1<i<T}\|\alpha_i\|_2+ \lambda*\|E\|_1, s.t. \|D_i\|_2= 1, 
%     (sequentially rank-1 fitting. The residual of the first column is always A.) 
%     [A, D, alpha] = LowRankDictionaryShrinkage_seq(A_ini, D, alpha, 1/mu); 



    
    

    Iteration = Iteration + 1;
    
    Error = X - A - E;
    Z = (Z + Error)/rho;  


    mu = min(mu*rho, mu_bar);

        
    %% stop Criterion    
    stopCriterion = norm(Error, 'fro') / d_norm;
    if stopCriterion < tol || iter >= maxIter
        converged = true;
    end    
    

     disp(['#Iter ' num2str(Iteration) ' r(A) ' num2str(size(D,2))...
            ' |E|_0 ' num2str(length(find(abs(E)>0)))...
            ' stopCriterion ' num2str(stopCriterion)]);


end


