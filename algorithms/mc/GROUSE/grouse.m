%function [U,R,err_reg] = grouse(I,J,S,numr,numc,maxrank,step_size,maxCycles,Uinit)
function [U,R,err_reg] = grouse(values,Indicator,numr,numc,maxrank,step_size,maxCycles,Uinit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  GROUSE (Grassman Rank-One Update Subspace Estimation) matrix completion code 
%  by Ben Recht and Laura Balzano, February 2010.
%
%  Given a sampling of entries of a matrix X, try to construct matrices U
%  and R such that U is unitary and UR' approximates X.  This code 
%  implements a stochastic gradient descent on the set of subspaces.
%
%  Inputs:
%       (I,J,S) index the known entries across the entire data set X. So we
%       know that for all k, the true value of X(I(k),J(k)) = S(k)
%
%       numr = number of rows
%       numc = number of columns
%           NOTE: you should make sure that numr<numc.  Otherwise, use the
%           transpose of X
%       
%       max_rank = your guess for the rank
%
%       step_size = the constant for stochastic gradient descent step size
%
%       maxCycles = number of passes over the data
%
%       Uinit = an initial guess for the column space U (optional)
%
%   Outputs:
%       U and R such that UR' approximates X.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matlab specific data pre-processing
%

% Form some sparse matrices for easier matlab indexing
%values = sparse(I,J,S,numr,numc);
%Indicator = sparse(I,J,1,numr,numc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Main Algorithm
%

if (nargin<9)
    % initialize U to a random r-dimensional subspace 
    U = orth(randn(numr,maxrank)); 
else
    U = Uinit;
end

err_reg = zeros(maxCycles*numc,1);


for outiter = 1:maxCycles,
    
    fprintf('Pass %d...\n',outiter);
    
    % create a random ordering of the columns for the current pass over the
    % data.
    col_order = randperm(numc);
    
for k=1:numc,
      
    % Pull out the relevant indices and revealed entries for this column
    idx = find(Indicator(:,col_order(k)));
    v_Omega = values(idx,col_order(k));
    U_Omega = U(idx,:);    

    
    % Predict the best approximation of v_Omega by u_Omega.  
    % That is, find weights to minimize ||U_Omega*weights-v_Omega||^2
    
    weights = U_Omega\v_Omega;
    norm_weights = norm(weights);
    
    % Compute the residual not predicted by the current estmate of U.

    residual = v_Omega - U_Omega*weights;       
    norm_residual = norm(residual);
    
    % This step-size rule is given by combining Edelman's geodesic
    % projection algorithm with a diminishing step-size rule from SGD.  A
    % different step size rule could suffice here...        
    
    sG = norm_residual*norm_weights;
    err_reg((outiter-1)*numc + k) = norm_residual/norm(v_Omega);
    t = step_size*sG/( (outiter-1)*numc + k );
   
    % Take the gradient step.    
    if t<pi/2, % drop big steps        
        alpha = (cos(t)-1)/norm_weights^2;
        beta = sin(t)/sG;

        step = U*(alpha*weights);
        step(idx) = step(idx) + beta*residual;

        U = U + step*weights';
    end 
end

end

% Once we have settled on our column space, a single pass over the data
% suffices to compute the weights associated with each column.  You only
% need to compute these weights if you want to make predictions about these
% columns.
fprintf('Find column weights...');
R = zeros(numc,maxrank);
for k=1:numc,     
    % Pull out the relevant indices and revealed entries for this column
    idx = find(Indicator(:,k));
    v_Omega = values(idx,k);
    U_Omega = U(idx,:);
    % solve a simple least squares problem to populate R
    R(k,:) = (U_Omega\v_Omega)';
end

