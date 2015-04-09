function [X,S,out] = nsa_v2_original(D,stdev,tol,optimal_X,optimal_S)

% IMPORTANT:
% PROPACK package is required. You may download from
% http://soi.stanford.edu/~rmunk/PROPACK/


%************************************************************
% stop_type: stopping criterion
%            1 - NSA stops when
%               |(X,S)-(Xp,Sp)|_F/(|(Xp,Sp)|_F+1)< tol*stdev
%            2 - NSA stops when
%               |X-Z|_F/|D|_F < tol
%************************************************************
% Written by Necdet Serhat Aybat, Columbia University
%
% 2)last modified on 27 October 2011
%       - a bug is fixed related to the first stopping criterion    
% 1)first created on 28 July 2011.
%************************************************************
    stop_type = 1;
    if nargin == 3
        output_flag = 0;
    elseif nargin == 5
        output_flag = 1;
    end
    addpath(genpath('./'));
    
    tic
    [n1,n2] = size(D);
    n_max = max(n1,n2);
    n_min = min(n1,n2);
    
    deltabar = sqrt(n_min+sqrt(8*n_min))*stdev;
    X=0;
    norm_two = lansvd(D, 1, 'L');
    mu = 1.25/norm_two;     % this one can be tuned
    mu_bar = mu * 1e8;
    coef = 1.65;            % this one can be tuned
     
    iter = 1; 
    flag=true;
    sv_array = [];
    %sv = round(n_min/10);
    sv = 10;
    sv_array = [sv_array, sv];
    
    lambda = lambda_search(abs(D(:)), mu, deltabar, 1/sqrt(n_max));
    rho = (mu+lambda)/(sqrt(n_max)*mu*lambda);
        
    S = sign(D).*max(abs(D)-rho,0);
    Z = lambda/(lambda+mu)*(D-S);
    Y = -mu*Z;
    
    while(flag==true)
        
        if choosvd(n_min, sv) == 1
            [U,Sigma,V] = lansvd(Z-Y/mu, sv, 'L');
        else
            [U,Sigma,V] = svd(Z-Y/mu, 'econ');
        end
        diagSigma = diag(Sigma);
        svp = length(find(diagSigma > 1/mu));
        if svp < sv
            sv = min(svp + 1, n_min);
            disp('svp<sv')
        elseif svp==sv
            sv = min(svp + round(0.04*n_min), n_min);
            disp('svp=sv')
        end
        sv_array = [sv_array, sv];
        X_old = X;
        X = U(:, 1:svp) * diag(diagSigma(1:svp) - 1/mu) * V(:, 1:svp)';
        
        %An alternative memory efficient way to compute X
        %X = bsxfun(@times,U(:, 1:svp),(diagSigma(1:svp) - 1/mu)') * V(:, 1:svp)';
        
        S_old = S;
        S = D - X - Y/mu;
        lambda = lambda_search(abs(S(:)), mu, deltabar, 1/sqrt(n_max));
        rho = (mu+lambda)/(sqrt(n_max)*mu*lambda);
        
        S = sign(S).*max(abs(S)-rho,0);
        Z = lambda/(lambda+mu)*(D-S)+ mu/(lambda+mu)*(X+Y/mu);
        
        disp(['Iteration no: ',num2str(iter),' --> |X-Z|_F/|D|_F=', num2str(norm(X-Z,'fro')/norm(D,'fro')),', |X+S-D|_F/|D|_F=', num2str(norm(X+S-D,'fro')/norm(D,'fro'))])
        if stop_type == 1
            stop_crit = sqrt(norm(X-X_old,'fro')^2+norm(S-S_old,'fro')^2)/(sqrt(norm(X_old,'fro')^2+norm(S_old,'fro')^2)+1);
            converged = stop_crit <= tol*stdev;
        elseif stop_type == 2
            stop_crit = norm(X-Z,'fro')/norm(D,'fro');
            converged = stop_crit <= tol;
        end
        if converged
            total_time = toc
            disp([num2str(iter),' itererations done'])
            flag=false;
        end
        Y = Y + mu*(X-Z);
        mu = min(mu*coef, mu_bar);
        iter=iter+1;
    end
    
    disp([num2str(iter-1),' total partial svd number'])
    treshold =1e-5;
    Sigma_X = svd(X,'econ');
    o1 = sum(Sigma_X);
    o2 = norm(S(:),1);
    o3 = norm(X+S-D,'fro');

    if output_flag == 0
        out.rel_err_X = norm(X-D,'fro')/norm(D,'fro');
        out.rel_err_S = norm(S-D,'fro')/norm(D,'fro');
        out.rel_err_noise = o3/norm(D,'fro');
        out.rank = sum(Sigma_X>3*stdev);
        
        display(' ')
        display('**********************************************************')
        display(['Nuclear Norm: ', num2str(o1), '  L1 Norm: ', num2str(o2)])
        display(['|X+S-D|_F: ', num2str(o3)])
        display(['|D-X|/|D|: ',num2str(out.rel_err_X)])
        display(['|D-S|/|D|: ',num2str(out.rel_err_S)])
        display(['Rank: ',num2str(out.rank)])
        display(['|X+S-D|/|D|: ',num2str(out.rel_err_noise)])
    else       
        Sigma_opt = svd(optimal_X, 'econ');
        nuc_opt = sum(Sigma_opt);
        l1_opt = norm(optimal_S(:),1);
        out.rel_err_X = norm(X-optimal_X,'fro')/norm(optimal_X,'fro');
        out.rel_err_S = norm(S-optimal_S,'fro')/norm(optimal_S,'fro');
        nuc_norm_rel_err = abs(nuc_opt-o1)/nuc_opt;
        L1_rel_err = abs(l1_opt-o2)/l1_opt;
        out.rel_err_noise = o3/norm(D,'fro');
        out.rank = sum(Sigma_X>3*stdev);

        display(' ')
        display('**********************************************************')
        display(['Nuclear Norm(X*): ', num2str(nuc_opt),'  Nuclear Norm: ', num2str(o1)])
        display(['L1 Norm(S*): ', num2str(l1_opt),'  L1 Norm: ', num2str(o2)])
        display(['  |X+S-D|_F: ', num2str(o3)])
        display(['Relative Error X: ',num2str(out.rel_err_X)])
        display(['Relative Error s: ',num2str(out.rel_err_S)])
        
        display(' ')
        ind_pos_X = find(Sigma_opt>treshold);
        ind_0_X = find(Sigma_opt<=treshold);
        inf_err_pos_X = norm(Sigma_X(ind_pos_X)-Sigma_opt(ind_pos_X),inf);
        inf_err_0_X = norm(Sigma_X(ind_0_X)-Sigma_opt(ind_0_X),inf);
        display(['Relative Error Nuclear Norm: ',num2str(nuc_norm_rel_err)])
        display(['max{|sigmaX(i)-sigmaX*(i)|: i s.t. sigmaX*(i)>0}: ',num2str(inf_err_pos_X)])
        display(['max{|sigmaX(i)-sigmaX*(i)|: i s.t. sigmaX*(i)=0}: ',num2str(inf_err_0_X)])
        
        display(' ')
        ind_pos_S = find(abs(optimal_S)>treshold);
        ind_0_S = find(abs(optimal_S)<=treshold);
        inf_err_pos_S = norm(S(ind_pos_S)-optimal_S(ind_pos_S),inf);
        inf_err_0_S = norm(S(ind_0_S)-optimal_S(ind_0_S),inf);
        display(['Relative Error L1 Norm: ',num2str(L1_rel_err)])
        display(['max{|S(ij)-S*(ij)|: ij s.t. S*(ij)>0}: ',num2str(inf_err_pos_S)])
        display(['max{|S(ij)-S*(ij)|: ij s.t. S*(ij)=0}: ',num2str(inf_err_0_S)])
        
        display(' ')
        display(['Rank: ',num2str(out.rank)])
        display(['|X+S-D|/|D|: ',num2str(out.rel_err_noise)])
    end
end