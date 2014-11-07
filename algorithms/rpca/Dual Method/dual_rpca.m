function [A_dual E_dual Y iter] = dual_rpca(D, lambda, tol, maxIter, LineSearchFlag, outputFile)

% This matlab code implements the dual method for RPCA.
%
% The Primal Robust PCA relaxation
% min \tau ( |A|_* + \lambda |E|_1 ) + 1/2 |(A,E)|_F^2
% subj  A+E = D
%
% The Dual problem
% max trace(D' * Y)
% subj max( |Y|_2, 1 \ \lambda |Y|_inf) <= 1
%
% References
% "Robust PCA: Exact Recovery of Corrupted Low-Rank Matrices via Convex Optimization", J. Wright et al., preprint 2009.
% "Fast Convex Optimization Algorithms for Exact Recovery of a Corrupted Low-Rank Matrix", Z. Lin et al.,  preprint 2009.
%
% Minming Chen, June 2009. Questions? v-minmch@microsoft.com ; 
% Arvind Ganesh (abalasu2@illinois.edu)
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing

addpath PROPACK;

[m n] = size(D);

if nargin < 2
    lambda = 1 / sqrt(m);
end

if nargin < 3
    tol = 2e-5 * norm( D, 'fro' );
elseif tol == -1
    tol = 2e-5 * norm( D, 'fro' );
end

if nargin < 4
    maxIter = 1000;
elseif maxIter == -1
    maxIter = 1000;
end

if nargin < 5
    LineSearchFlag = 0;
elseif LineSearchFlag == -1
    LineSearchFlag = 0;
end

if nargin > 5
    fid = fopen(outputFile,'a+t') ;
end

% initialize
Y = sign(D);
norm_two = lansvd_d1(Y);
norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;
norm_two = norm_two / dual_norm;
norm_inf = norm_inf / dual_norm;
obj_v = D(:)' * Y(:);

% projection
A_dual = zeros( m, n);
E_dual = zeros( m, n);
eps = 1e-4;
tolProj = 1e-8 * norm( D, 'fro' );
epsProjection = 1e-1;
t = 1;

% linesearch
delta = 0.1;
K = 7;
memoStepsize = ones(1,K) * 0.1;
beta = 0.6;

iter = 0;
converged = false;
while ~converged       
    %% compute Z : the projction of D onto the normal cone of Y
    % get the search direction D - Z
    
    iter = iter + 1;
    if norm_two < norm_inf - eps && iter < maxIter
        
        threshold = norm(Y(:), inf) * ( 1 - epsProjection );
        Z = max( D .* (Y > threshold), 0) + min( D .* (Y < -threshold), 0);
        
    else
        
        t = max( round(t * 1.1), t + 1);
        if choosvd( n, t) == 1
            [u s v] = lansvd( Y, t, 'L');
        else
            [u s v] = svd( Y, 'econ');
        end
        ds = diag(s);
        t = max( find( ds >= ds(1) * (1 - 1e-2) ) );
        
        if norm_two > norm_inf + eps && iter < maxIter
            
            D_bar = u(:, 1:t)' * D * v(:, 1:t);
            [S J] = eig(( D_bar + D_bar') / 2);
            temp = S * max( J, 0) * S';
            Z = u(:, 1:t) * temp * v(:, 1:t)';
            
        else
            
            convergedProjection = false;
            A_dual = zeros( m, n);
            E_dual = zeros( m, n);
            proj = 0;    
            threshold = norm(Y(:), inf) * ( 1 - epsProjection );
            while ~convergedProjection
                Z = D - A_dual;
                Z = max( Z .* (Y > threshold), 0) + min( Z .* (Y < -threshold), 0);
                D_bar = u(:, 1:t)' * ( D - Z) * v(:, 1:t);
                [S J] = eig(( D_bar + D_bar') / 2);
                temp = S * max( J, 0) * S';
                X = u(:, 1:t) * temp * v(:, 1:t)';
    
                if norm( Z - E_dual, 'fro') < tolProj && norm( X - A_dual, 'fro') < tolProj
                    convergedProjection = true;
                    E_dual = Z;
                    A_dual = X;
                    Z = Z + X;
                else
                    E_dual = Z;
                    A_dual = X;
                end
                
                
                proj = proj + 1;    
                if proj > 50
                    convergedProjection = true;
                end
            end
            
        end
        
    end

    %% linesearch to find max trace(D'*(Y+delta*(D-Z)))/J(Y+delta*(D-Z))
    
    Z = D - Z;
    a = D(:)' * Y(:);
    b = D(:)' * Z(:);
    c = Z(:)' * Z(:);
    
    if LineSearchFlag == 0
        
        % nonexact linesearch
        stepsize = max( 1.3 * median(memoStepsize), 1e-4 );
        convergedLineSearchLike = 0;
        numTrialPoint = 0;
        while ~convergedLineSearchLike
            X = Y + stepsize * Z;
            norm_two = lansvd_d1(X);
            norm_inf = norm(X(:), inf) / lambda;
            dual_norm = max(norm_two, norm_inf);
            tempv = (a + b * stepsize) / dual_norm;
            diff = tempv - a - stepsize / 2 * c;
            if diff > 0 || numTrialPoint >= 50
                convergedLineSearchLike = 1;
                obj_v = tempv;
                norm_two = norm_two / dual_norm;
                norm_inf = norm_inf / dual_norm;
                Y = X / dual_norm;
                delta = stepsize;
            else
                stepsize = stepsize * beta;
            end
            numTrialPoint = numTrialPoint + 1;
        end        
        memoStepsize(1, mod(iter,K)+1) = delta; 
        
    elseif LineSearchFlag == 1
                
        epsLineSearch = 1e-10;
        stepsize = max(abs(median(memoStepsize)),1e-4);
    
        currentPosition = 0;
        direction = 1;
        currentValue = a;
        position_trace = 0;
        value_trace = currentValue;    
        numTrialPoint = 1;
        convergedLineSearch = 0;
        
        while ~convergedLineSearch
        
            nextPosition = currentPosition + direction * stepsize;
            findnp = find(abs(position_trace - nextPosition) < 1e-1 * epsLineSearch);
            if size(findnp , 2) > 0
                nextValue = value_trace(findnp);
            else
                X = Y + nextPosition * Z;
                norm_two = lansvd_d1(X);
                norm_inf = norm(X(:), inf) / lambda;
                dual_norm = max( norm_two, norm_inf );
                nextValue = ( a + nextPosition * b ) / dual_norm ;
                position_trace = [position_trace nextPosition];
                value_trace = [value_trace nextValue];
                numTrialPoint = numTrialPoint + 1;
            end
        
            if nextValue <= currentValue
                direction = -direction;
                stepsize = stepsize / 2;
            else
                currentPosition = nextPosition;
                currentValue = nextValue;
            end
        
            if stepsize < epsLineSearch || numTrialPoint >= 50 || stepsize < 1e-1 * currentPosition
                delta = currentPosition;
                obj_v = currentValue;
                Y = Y + delta * Z;
                norm_two = lansvd_d1(Y);
                norm_inf = norm(Y(:), inf) / lambda;
                dual_norm = max( norm_two, norm_inf);
                Y = Y / dual_norm;
                norm_two = norm_two / dual_norm;
                norm_inf = norm_inf / dual_norm;
                convergedLineSearch = 1;
            end        
        end 
        memoStepsize(1, mod(iter,K)+1) = delta;
        
    end
        
    %% stop Criterion    
    stopCriterion = norm( D - A_dual - E_dual, 'fro' );
    if stopCriterion < tol 
        converged = true;
    end    
    
    if mod( iter, 10) == 0
        disp([num2str(iter) ' Y_principal ' num2str(t)...
            ' |E|_0 ' num2str(length(find(abs(E_dual)>0)))...
            ' objvalue ' num2str(obj_v) ' |D-A-E|_F ' num2str(stopCriterion)]);
        if nargin > 5
            fprintf(fid, '%s\n', [num2str(iter) ' Y_principal ' num2str(t)... 
            ' |E|_0 ' num2str(length(find(abs(E_dual)>0)))...
            ' objvalue ' num2str(obj_v) ' |D-A-E|_F ' num2str(stopCriterion)]);
        end
    end    
    
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end

if nargin > 5
    fclose(fid);
end

