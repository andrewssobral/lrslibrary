function [A_hat,E_hat,numIter] = proximal_gradient_rpca(D, lambda,...
    maxIter, tol, lineSearchFlag, ...
    continuationFlag, eta, mu, outputFileName )

% D - m x n matrix of observations/data (required input)
% lambda - weight on sparse error term in the cost function (required input)
%
% tol - tolerance for stopping criterion.
%     - DEFAULT 1e-7 if omitted or -1.
% maxIter - maximum number of iterations
%         - DEFAULT 10000, if omitted or -1.
% lineSearchFlag - 1 if line search is to be done every iteration
%                - DEFAULT 0, if omitted or -1.
% continuationFlag - 1 if a continuation is to be done on the parameter mu
%                  - DEFAULT 1, if omitted or -1.
% eta - line search parameter, should be in (0,1)
%     - ignored if lineSearchFlag is 0.
%     - DEFAULT 0.9, if omitted or -1.
% mu - relaxation parameter
%    - ignored if continuationFlag is 1.
%    - DEFAULT 1e-3, if omitted or -1.
% outputFileName - Details of each iteration are dumped here, if provided.
%
% [A_hat, E_hat] - estimates for the low-rank part and error part, respectively
% numIter - number of iterations until convergence
%
%
% References
% "Robust PCA: Exact Recovery of Corrupted Low-Rank Matrices via Convex Optimization", J. Wright et al., preprint 2009.
% "An Accelerated Proximal Gradient Algorithm for Nuclear Norm Regularized Least Squares problems", K.-C. Toh and S. Yun, preprint 2009.
%
% Arvind Ganesh, Summer 2009. Questions? abalasu2@illinois.edu
%
% Copyright: Perception and Decision Laboratory
%			University of Illinois, Urbana-Champaign

if nargin < 2
    error('Too few arguments') ;
end

if nargin < 4
    tol = 1e-6 ;
elseif tol == -1
    tol = 1e-6 ;
end

if nargin < 3
    maxIter = 1000 ;
elseif maxIter == -1
    maxIter = 1000 ;
end


if nargin < 5
    lineSearchFlag = 0 ;
elseif lineSearchFlag == -1
    lineSearchFlag = 0 ;
end

if nargin < 6
    continuationFlag = 1 ;
elseif continuationFlag == -1 ;
    continuationFlag = 1 ;
end

if lineSearchFlag
    if nargin < 7
        eta = 0.9 ;
    elseif eta == -1
        eta = 0.9 ;
    end
    
    if ( nargin > 6 && ( eta <= 0 || eta >= 1 ) )
        disp('Line search parameter out of bounds, switching to default 0.9') ;
        eta = 0.9 ;
    end
end

if ~continuationFlag
    if nargin < 8
        mu = 1e-3 ;
    elseif mu == -1
        mu = 1e-3 ;
    end
end

if nargin > 8
    fid = fopen(outputFileName,'w') ;
end

DISPLAY_EVERY = 20 ;
DISPLAY = 1 ;

maxLineSearchIter = 200 ; % maximum number of iterations in line search per outer iteration

%% Initializing optimization variables

[m,n] = size(D) ;

t_k = 1 ; % t^k
t_km1 = 1 ; % t^{k-1}

tau_0 = 2 ; % square of Lipschitz constant for the RPCA problem

X_km1_A = zeros(m,n) ; X_km1_E = zeros(m,n) ; % X^{k-1} = (A^{k-1},E^{k-1})
X_k_A = zeros(m,n) ; X_k_E = zeros(m,n) ; % X^{k} = (A^{k},E^{k})


if continuationFlag
    if lineSearchFlag
        mu_0 = eta * norm(D) ;
    else
        mu_0 = norm(D) ;
    end
    
    mu_k = 0.99*mu_0 ;
    
    mu_bar = 1e-9 * mu_0 ;
else
    mu_k = mu ;
    
end

tau_k = tau_0 ;

converged = 0 ;
numIter = 0 ;

stoppingCriterionOld = -inf ;
stagnationEpsilon = 1e-6 ;

% tol = 1e-6 ;
oldCost = -inf ;
cost = [] ;
cost1 = [] ;

mu_path = mu_k ;

%% Start main loop

while ~converged
    
    Y_k_A = X_k_A + ((t_km1 - 1)/t_k)*(X_k_A-X_km1_A) ;
    Y_k_E = X_k_E + ((t_km1 - 1)/t_k)*(X_k_E-X_km1_E) ;
    
    if ~lineSearchFlag
        
        G_k_A = Y_k_A - (1/tau_k)*(Y_k_A+Y_k_E-D) ;
        G_k_E = Y_k_E - (1/tau_k)*(Y_k_A+Y_k_E-D) ;
        
        [U,S,V] = svd(G_k_A,'econ') ;
        diagS = diag(S) ;
        
        X_kp1_A = U * diag(pos(diagS-mu_k/tau_k)) * V';
        X_kp1_E = sign(G_k_E) .* pos( abs(G_k_E) - lambda*mu_k/tau_k );
        
        rankA  = sum(diagS>mu_k/tau_k);
        cardE = sum(sum(double(abs(X_kp1_E)>0)));
        
    else
        
        convergedLineSearch = 0 ;
        numLineSearchIter = 0 ;
        
        tau_hat = eta * tau_k ;
        
        while ~convergedLineSearch
            
            temp = (1/tau_hat)*(Y_k_A + Y_k_E - D) ;
            G_A = Y_k_A - temp ;
            G_E = Y_k_E - temp ;
            
            [U,S,V] = svd(G_A,'econ') ;
            diagS = diag(S) ;
            SG_A = U * diag(pos(diagS-mu_k/tau_hat)) * V';
            SG_E = sign(G_E) .* pos( abs(G_E) - lambda*mu_k/tau_hat );
            
            F_SG = 0.5*norm(D-SG_A-SG_E,'fro')^2 ;
            Q_SG_Y = 0.5*tau_hat*norm([SG_A,SG_E]-[G_A,G_E],'fro')^2 + ...
                (0.5-1/tau_hat)*norm(D-Y_k_A-Y_k_E,'fro')^2    ;
            
            if F_SG <= Q_SG_Y
                tau_k = tau_hat ;
                convergedLineSearch = 1 ;
            else
                tau_hat = min(tau_hat/eta,tau_0) ;
            end
            
            numLineSearchIter = numLineSearchIter + 1 ;
            
            if ~convergedLineSearch && numLineSearchIter >= maxLineSearchIter
                disp('Stuck in line search') ;
                convergedLineSearch = 1 ;
            end
            
        end
        
        X_kp1_A = SG_A ;
        X_kp1_E = SG_E ;
        
        rankA  = sum(diagS>mu_k/tau_hat);
        cardE = sum(sum(double(abs(X_kp1_E)>0)));
        
        
    end
    
    t_kp1 = 0.5*(1+sqrt(1+4*t_k*t_k)) ;
    
    temp = X_kp1_A + X_kp1_E - Y_k_A - Y_k_E ;
    S_kp1_A = tau_k*(Y_k_A-X_kp1_A) + temp ;
    S_kp1_E = tau_k*(Y_k_E-X_kp1_E) + temp ;
    
    
    stoppingCriterion = norm([S_kp1_A,S_kp1_E],'fro')/(tau_k*max(1,norm([X_kp1_A,X_kp1_E],'fro'))) ;
    
    if stoppingCriterion <= tol
        converged = 1 ;
    end
    
    if continuationFlag
        mu_k = max(0.9*mu_k,mu_bar) ;        
    end
    
    t_km1 = t_k ;
    t_k = t_kp1 ;
    X_km1_A = X_k_A ; X_km1_E = X_k_E ;
    X_k_A = X_kp1_A ; X_k_E = X_kp1_E ;
    
    numIter = numIter + 1 ;
    
    if DISPLAY && mod(numIter,DISPLAY_EVERY) == 0
        disp(['Iteration ' num2str(numIter) '  rank(A) ' num2str(rankA) ...
            ' ||E||_0 ' num2str(cardE)]) ;         
    end
    
    if nargin > 8
        fprintf(fid, '%s\n', ['Iteration ' num2str(numIter) '  rank(A)  ' num2str(rankA) ...
            '  ||E||_0  ' num2str(cardE) '  Stopping Criterion   ' ...
            num2str(stoppingCriterion)]) ;
    end
    
    
    if ~converged && numIter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;
    end
    
end

A_hat = X_k_A ;
E_hat = X_k_E ;

if nargin > 8
    fprintf(fid,'\n\n') ;
    if lineSearchFlag
        fprintf(fid,'Line search ON\n') ;
    else
        fprintf(fid,'Line search OFF\n') ;
    end
    if continuationFlag
        fprintf(fid,'Continuation ON\n') ;
    else
        fprintf(fid,'Continuation OFF\n') ;
    end
    fclose(fid) ;
end