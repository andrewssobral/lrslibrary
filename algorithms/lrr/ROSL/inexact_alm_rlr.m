function [alpha, E, A] = inexact_alm_rlr(X, K, D, tol, maxIter) 

% This matlab code implements the inexact augmented Lagrange multiplier 
% method for Robust linear regresion with known alpha
%
% alpha- k x n matrix of known coefficients
%
% X - m x n matrix of observations/data (required input)
%
% tol - tolerance for stopping criterion.
%
% maxIter - maximum number of iterations
%  
% 
% min_D \|E\|_1, s.t. X = D*alpha +E;
% 
% The laplacian function: L(A,E,Y,u) =  |E|_1  + \mu/2* |X- D*alpha-E+Y/mu|_F^2;
%
%
% Xianbiao Shu (xshu2@illinois.edu)
% Copyright: Mitsubishi Electric Research Lab
% Reference: X. Shu, F. Porikli, N. Ahujia "Robust Orthonormal Subspace Learning: Efficient Recovery of Corrupted Low-rank Matrices". CVPR 2014; 

addpath PROPACK;


[m, n] = size(X);

% initialize
% alpha   % Low-rank coeffs

E = zeros(size(X)); % sparse innovation
A = X;


Z = 0; % Z = Y/mu

norm_inf = norm(X(:), inf); 
d_norm = norm(X, 'fro');

mu = 10*5e-2/norm_inf; % this parameter can be tuned
rho = 1.5;       %1.3    % this parameter can be tuned
mu_bar = mu * 1e7;



iter = 0;
Iteration = 0;
converged = false;


while ~converged    

    iter = iter + 1;
    
    E_temp = X - A + Z; 

    


    %% Intensity threshold
    E = max(abs(E_temp) - 1/mu, 0).*sign(E_temp);   




    
    %% Sparse Representation
    A_ini = X -E + Z;
    
    

%%   given D and A
    
    [U, S, V] = svd(D, 0); 
    diag_S = diag(S);
    SV = length(diag_S(diag_S>0));
    
    alpha = V(:, 1:SV)*diag(1./diag_S(1:SV))*(U(:, 1:SV))'*A_ini;
    A = D*alpha;
    
    
        

    
    Iteration = Iteration + 1;
    
    Error = X - A - E;
    Z = (Z + Error)/rho;  


    mu = min(mu*rho, mu_bar);

        
    %% stop Criterion    
    stopCriterion = norm(Error, 'fro') / d_norm;
    if stopCriterion < tol || iter >= maxIter
        converged = true;
    end    
    

     disp(['#Iter ' num2str(Iteration)  ...
            ' |E|_0 ' num2str(length(find(abs(E)>0)))...
            ' stopCriterion ' num2str(stopCriterion)]);


end


