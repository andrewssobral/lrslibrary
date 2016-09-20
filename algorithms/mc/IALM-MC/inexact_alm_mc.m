function [A iter svp] = inexact_alm_mc(D, tol, maxIter, rho)

% Oct 2009
% This matlab code implements the inexact augmented Lagrange multiplier 
% method for Matrix Completion.
%
% D - m x n matrix of observations/data (required input)
%
% tol - tolerance for stopping criterion.
%     - DEFAULT 1e-7 if omitted or -1.
%
% maxIter - maximum number of iterations
%         - DEFAULT 1000, if omitted or -1.
%
% Model: 
%     min |A|_*
%     subj A + E = D, \pi_{\Omega}(E) = 0
%
% Algorithm:
%
% Initialize A,E,Y,u
% while ~converged 
%   minimize (inexactly, update A and E only once)
%     L(A,E,Y,u) = |A|_* + <Y,D-A-E> + mu/2 * |D-A-E|_F^2;
%   Y = Y + \mu * (D - A - E);
%   update \mu;
% end
%
% Minming Chen (cmmfir@gmail.com);
% Zhouchen Lin (zhoulin@microsoft.com; zhouchenlin@gmail.com)
%
% Reference: Zhouchen Lin, Minming Chen, and Yi Ma, The Augmented Lagrange Multiplier Method 
%for Exact Recovery of Corrupted Low-Rank Matrix, http://perception.csl.illinois.edu/matrix-rank/Files/Lin09-MP.pdf
%
% Copyright: Microsoft Research Asia, Beijing

%clear global;
global A Sparse_Z;

%addpath PROPACK;

[m,n] = size(D);

if nargin < 2
    tol = 1e-4;
elseif tol == -1
    tol = 1e-4;
end

if nargin < 3
    maxIter = 1000;
elseif maxIter == -1
    maxIter = 1000;
end

% read sparse matrix D
[m n] = size(D);
[I J V] = find(D);
p = length(I);
col = [0; find(diff(J)); p];
Sparse_Z = sparse(D);
clear D;

% initialize
sv = 5;
svp = sv;
A.U = zeros(m, sv);
A.V = zeros(n, sv);
Y = zeros(p, 1);
Z = zeros(p, 1); % Z indicates the nonzeros of \pi_{\Omega}(A)
norm_D = norm(V, 'fro');
rho_s = p / (m * n);
mu = 0.3/(lansvd('Axz','Atxz',m,n,1,'L'));
mu_bar = mu * 1e20;
tol2 = (1e-5) * sqrt(1 - rho_s) / sqrt(rho_s);
if nargin < 4
    rho = 1.8588*rho_s + 1.2172;
elseif rho == -1
    rho = 1.8588*rho_s + 1.2172;
end

% Iterations 
iter = 0;
converged1 = false;
converged2 = false;
stopCriterion1 = 1;
stopCriterion2 = 1;
while ~converged1 || ~converged2 
    iter = iter + 1;    
    converged1 = false;
    converged2 = false;
    
    A_old = A;
    Z_old = Z;
    
    %% compute (partial) SVD for update A
    Sparse_Z = spconvert([I,J,V-Z+1/mu*Y; m,n,0]);
    if stopCriterion1 > 10 * tol
        options.tol = 10*tol;
    else
        options.tol = min(0.1*tol, 0.01/mu);
    end
    [A.U,S,A.V] = lansvd('Axz','Atxz',m,n,sv,'L',options);
        
    %% get the rank of A.
    diagS = diag(S);
    diagS = diagS(1:sv);
    svn = length(find(diagS > 1/mu));
    svp = svn;
    
    %% truncation
    ratio = diagS(1:end-1)./diagS(2:end);
    [max_ratio, max_idx] = max(ratio);
    if max_ratio > 2
        svp = min(svn, max_idx);
    end
    if svp < sv 
        sv = min(svp + 1, min(m,n));
    else
        sv = min(svp + 10, min(m,n));
    end
        
    %% update A
    sqrtds = sqrt(diagS(1:svp) - 1/mu);
    A.U = A.U(:, 1:svp) * diag(sqrtds);
    A.V = A.V(:, 1:svp) * diag(sqrtds);
    
    %% update Z and Y
    Z = UVtOmega(A.U,A.V,I,J,col);
    Y = Y + mu * (V - Z);
        
    %% The first stop criterion is ||D-A-E|| / ||D|| < tol;
    %% Note that \pi_{\Omega}(E) = 0 and \pi_{\bar{\Omega}}(E) = \pi_{\Omega}(A)
    stopCriterion1 = norm(V - Z, 'fro') / norm_D;
    if stopCriterion1 < tol
        converged1 = true;
    end    
    
    %%The second stop criterion is mu * ||E_{k+1} - E_{k}|| / ||D|| < tol2
    norm_diffA = sum(sum((A_old.U'*A_old.U).*(A_old.V'*A_old.V))) +...
        sum(sum((A.U'*A.U).*(A.V'*A.V))) - 2*sum(sum((A.U'*A_old.U).*(A.V'*A_old.V)));
    stopCriterion2 = min(mu,sqrt(mu)) * sqrt(norm_diffA - norm(Z - Z_old, 'fro')^2) / norm_D; 
    if stopCriterion2 < tol2
        converged2 = true;
    end

    %% update mu   
    if converged2 
        mu = min(rho*mu, mu_bar);
    end

    %% display
    if mod( iter, 1) == 0
        disp(['#svd ' num2str(iter) ' r(A) ' num2str(svp)...
            ' svn ' num2str(svn) ' max_idx ' num2str(max_idx) ' sv ' num2str(sv) ...
            ' stopCriterion1 ' num2str(stopCriterion1)...
            ' stopCriterion2 ' num2str(stopCriterion2)...
            ' mu ' num2str(mu)]);
    end    
    
    %% Maximum iterations reached
    if (~converged1 || ~converged2) && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged1 = 1;       
        converged2 = 1;       
    end   
end

