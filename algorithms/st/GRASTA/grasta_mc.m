function [ U_hat, R,Outliers ] = grasta_mc( I,J,S, numr,numc,maxCycles,CONVERGE_LEVLE,OPTIONS)
%  GRASTA (Grassmannian Robust Adaptive Subspace Tracking Algorithm) robust 
%  matrix completion code   
%  by Jun He and Laura Balzano, Sept. 2011.
%
%   Online Robust Subspace Tracking from Partial Information
%       http://arxiv.org/abs/1109.3827
%
% Inputs: 
%
%       (I,J,S) index the known entries across the entire data set X. So we
%       know that for all k, the true value of X(I(k),J(k)) = S(k)
%
%       numr = number of rows
%       numc = number of columns
%           NOTE: you should make sure that numr<numc.  Otherwise, use the
%           transpose of X
%
% Outputs: 
%       U and R such that UR' approximates X.
%
%
% Form some sparse matrices for easier matlab indexing
values    = sparse(I,J,S,numr,numc);
Indicator = sparse(I,J,1,numr,numc);

status = struct();  % maintain GRASTA running status
status.init  = 0;   % will be set 1 once GRASTA start working

OPTS  = struct(); % initial a empty struct for OPTS
U_hat = zeros(1); % U_hat will be initialized in GRASTA


for outIter = 1 : maxCycles,
    
    col_order = randperm(numc);
    for k = 1 : numc,
        idx = find(Indicator(:,col_order(k)));
        v_Omega = values(idx,col_order(k));
        
        if length(idx) < OPTIONS.RANK * 1,
            continue;
        end
                
        [U_hat, status, OPTS] = grasta_stream(v_Omega, idx, U_hat, status, OPTIONS, OPTS);
                
    end

    if status.level >= CONVERGE_LEVLE,
%         fprintf('Pass %d/%d, reach the convergence level - %d...\n',outIter, maxCycles,status.level);
        break;
    end
    
%     fprintf('Pass %d/%d ......\n',outIter, maxCycles);
end


% OPTS2 used for recovering R
OPTS2 = OPTS;

R = zeros(numc,OPTIONS.RANK);
Outliers =zeros(numc,numr);
for k=1:numc,
    idx = find(Indicator(:,k));
    v_Omega = values(idx,k);
    
    if length(idx) < OPTIONS.RANK * 1,
        continue;
    end
        
    U_Omega = U_hat(idx,:);
    
    if OPTIONS.USE_MEX,
        [s, w, ~] = mex_srp(U_Omega, v_Omega, OPTS2);
    else
        [s, w, ~] = admm_srp(U_Omega, v_Omega, OPTS2);
    end

    R(k,:) = w';
    Outliers(k,idx) = s;
end

end
