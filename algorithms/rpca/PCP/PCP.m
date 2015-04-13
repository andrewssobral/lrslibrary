function [L,S] = PCP(M,lam,tol)
%
% [L,S] = PCP(M,lam,opts)

% This code solves the following model
%
% min_A { lam*||S(:)||_1 + ||L||_* }
% s.t. M = S+L

% where M is the data matrix, which will be decomposed into
% S sparse matrix S and low-rank matrix L.
%
% lam -- S small positive parameter

%% parameter setting
beta = .25/mean(abs(M(:))); % 
maxit = 1000;

%% initialization
[m,n] = size(M);
S = zeros(m,n);
L = zeros(m,n);
Lambda = zeros(m,n); % the dual variable

% main
for iter = 1:maxit
    
    nrmLS = norm([S,L],'fro');
    % dS, dL record the change of S and L, only used for stopping criterion
    
    %% S - subproblem 
    % S = argmin_A  lam*||S||_1 - <Lambda, S+L-M> + (beta/2) * ||S+L-M||.^2
    % Define element wise softshinkage operator as 
    %     softshrink(z; gamma) = sign(z).* max(abs(z)-gamma, 0);
    % S has closed form solution: S=softshrink(Lambda/beta + M - L; lam/beta)
    % (see my slide page 42 Equation (66).
    X = Lambda / beta + M;
    Y = X - L;
    dS = S;
    S = sign(Y) .* max(abs(Y) - lam/beta, 0); % softshinkage operator
    dS = S - dS;
    
    %% L - subproblem
    % L = argmin_B ||L||_* -<Lambda, S+L-M> + + (beta/2) * ||S+L-M||.^2
    % L has closed form solution (singular value thresholding)
    % see my slide page 42, Equation (65).
    Y = X - S;
    dL = L;
 
    %[U,D,V] = svd(Y,'econ'); % use 'econ' is more efficient especially when Y is large
    [U,D,V] = svdecon(Y); % fastest
    
    VT=V';
    D = diag(D);
    ind = find(D > 1/beta);
    D = diag(D(ind) - 1/beta);
    L = U(:,ind) * D * VT(ind,:);
    dL = L - dL;
    
    %% stopping criterion
    RelChg = norm([dS,dL],'fro') / (1 + nrmLS);
    %fprintf('Iter %d, RelChg %4.2e \n',iter,RelChg);
    if RelChg < tol, break; end
    
    %% Update Lambda (dual variable)
    Lambda = Lambda - beta * (S + L - M);
end
end
