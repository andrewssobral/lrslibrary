function [X,S,out] = nsa_original(D,stdev,tol,denoise_flag,optimal_X,optimal_S)

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
% Written by Necdet Serhat Aybat, Penn State University
%
% 3)denoising added as an option on 15 October 2012
% 2)last modified on 27 October 2011
%       - a bug is fixed related to the first stopping criterion    
% 1)first created on 28 July 2011.
%************************************************************
    stop_type = 1;
    if nargin <= 4
        output_flag = 0;
    elseif nargin == 6
        output_flag = 1;
    end
    addpath(genpath('./'));
    display('**********************************************************')
    display('LEGEND')
    display('**********************************************************')
    display('(X,S): iterates computed at current iteration')
    display('(Xp,Sp): iterates computed at previous iteration')
    if nargin == 6
        display('(X*,S*): Low-rank and sparse components of D=X*+S*+E*')
        display('relX=|X-X*|_F/|X*|_F, relS=|S-S*|_F/|S*|_F')
    end
    display('svd: total number of partial svds computed')
    display('lsv: total number of leading singular values computed / svd')
    display('cpu: total run time')
    display('**********************************************************')
    
    tic
    [n1,n2] = size(D);
    n_max = max(n1,n2);
    n_min = min(n1,n2);
    
    delta = sqrt(n_min+sqrt(8*n_min))*stdev;
    X=0;
    norm_two = lansvd(D, 1, 'L');
    rho = 1.25/norm_two;     % this one can be tuned
    rho_bar = rho * 1e8;
    kappa = 1.65;            % this one can be tuned
     
    iter = 1; 
    flag=true;
    sv_array = [];
    sv = round(n_min/10);
       
    theta = theta_search(abs(D(:)), rho, delta, 1/sqrt(n_max));
    c = (rho+theta)/(sqrt(n_max)*rho*theta);
        
    S = sign(D).*max(abs(D)-c,0);
    Z = theta/(theta+rho)*(D-S);
    Y = -rho*Z;
    
    while(flag==true)
        
        if choosvd(n_min, sv) == 1
            [U,Sigma,V] = lansvd(Z-Y/rho, sv, 'L');
        else
            [U,Sigma,V] = svd(Z-Y/rho, 'econ');
	     sv = n_min;
        end
        sv_array = [sv_array, sv];

        Sigma = diag(Sigma);
        svp = length(find(Sigma > 1/rho));
        if svp < sv
            sv = min(svp + 1, n_min);
        elseif svp==sv
            sv = min(svp + round(0.04*n_min), n_min);
        end
        
        X_old = X;
        X = U(:, 1:svp) * diag(Sigma(1:svp) - 1/rho) * V(:, 1:svp)';
        
        %An alternative memory efficient way to compute X
        %X = bsxfun(@times,U(:, 1:svp),(Sigma(1:svp) - 1/rho)') * V(:, 1:svp)';
        
        S_old = S;
        S = D - X - Y/rho;
        theta = theta_search(abs(S(:)), rho, delta, 1/sqrt(n_max));
        c = (rho+theta)/(sqrt(n_max)*rho*theta);
        
        S = sign(S).*max(abs(S)-c,0);
        Z = theta/(theta+rho)*(D-S)+ rho/(theta+rho)*(X+Y/rho);
        
        if stop_type == 1
            stop_crit = sqrt(norm(X-X_old,'fro')^2+norm(S-S_old,'fro')^2)/(sqrt(norm(X_old,'fro')^2+norm(S_old,'fro')^2)+1);
            disp(['Iteration no: ',num2str(iter),' --> |(X,S)-(Xp,Sp)|_F/(1+|(Xp,Sp)|_F)=', num2str(stop_crit)])
            converged = stop_crit <= tol*stdev;
        elseif stop_type == 2
            stop_crit = norm(X-Z,'fro')/norm(D,'fro');
            disp(['Iteration no: ',num2str(iter),' --> |X-Z|_F/|D|_F=', num2str(stop_crit)])
            converged = stop_crit <= tol;
        end
        if converged
            total_time = toc;
            flag=false;
        end
        Y = Y + rho*(X-Z);
        rho = min(rho*kappa, rho_bar+iter);
        iter=iter+1;
    end
    
    if denoise_flag == 1
        [X,S] = post_noise_removal(X,D,stdev);
    end
    noise_treshold =1e-5;
    nrm_D=norm(D,'fro');
    Sigma_X = svd(X,'econ');
    o1 = sum(Sigma_X);
    o2 = norm(S(:),1);
    o3 = norm(X+S-D,'fro');
    out.rel_err_noise = o3/nrm_D;
    out.rank = sum(Sigma_X>3*stdev);
    if output_flag == 0 
        display('**********************************************************')
        display('NSA Report:')
        display(['svd = ',num2str(iter-1)])
        display(['lsv = ',num2str(mean(sv_array))])
        display(['cpu = ',num2str(total_time)])
        display(['Rank: ',num2str(out.rank)])
        display(['|X+S-D|_F/|D|_F: ',num2str(out.rel_err_noise)])
    else       
        Sigma_opt = svd(optimal_X, 'econ');
        nuc_opt = sum(Sigma_opt);
        l1_opt = norm(optimal_S(:),1);
        out.rel_err_X = norm(X-optimal_X,'fro')/norm(optimal_X,'fro');
        out.rel_err_S = norm(S-optimal_S,'fro')/norm(optimal_S,'fro');
        nuc_norm_rel_err = abs(nuc_opt-o1)/nuc_opt;
        L1_rel_err = abs(l1_opt-o2)/l1_opt;
        display('**********************************************************')
        display('NSA Report:')
        display(['svd = ',num2str(iter-1)])
        display(['lsv = ',num2str(mean(sv_array))])
        display(['cpu = ',num2str(total_time)])
        display(['relX: ',num2str(out.rel_err_X)])
        display(['relS: ',num2str(out.rel_err_S)])
        display(['Rank: ',num2str(out.rank)])
        display(['|X+S-D|_F/|D|_F: ',num2str(out.rel_err_noise)])
        display('**********************************************************')
        display('Additional Statistics')
        display('**********************************************************')
        display(['Nuclear Norm(X*): ', num2str(nuc_opt),'  Nuclear Norm(X): ', num2str(o1)])
        display(['L1 Norm(S*): ', num2str(l1_opt),'  L1 Norm(S): ', num2str(o2)])
        
        display(' ')
        ind_pos_X = find(Sigma_opt>noise_treshold);
        ind_0_X = find(Sigma_opt<=noise_treshold);
        inf_err_pos_X = norm(Sigma_X(ind_pos_X)-Sigma_opt(ind_pos_X),inf);
        inf_err_0_X = norm(Sigma_X(ind_0_X)-Sigma_opt(ind_0_X),inf);
        display(['Relative Error Nuclear Norm: ',num2str(nuc_norm_rel_err)])
        display(['max{|sigmaX(i)-sigmaX*(i)|: i s.t. sigmaX*(i)>0}: ',num2str(inf_err_pos_X)])
        display(['max{|sigmaX(i)-sigmaX*(i)|: i s.t. sigmaX*(i)=0}: ',num2str(inf_err_0_X)])
        
        display(' ')
        ind_pos_S = find(abs(optimal_S)>noise_treshold);
        ind_0_S = find(abs(optimal_S)<=noise_treshold);
        inf_err_pos_S = norm(S(ind_pos_S)-optimal_S(ind_pos_S),inf);
        inf_err_0_S = norm(S(ind_0_S)-optimal_S(ind_0_S),inf);
        display(['Relative Error L1 Norm: ',num2str(L1_rel_err)])
        display(['max{|S(ij)-S*(ij)|: ij s.t. S*(ij)>0}: ',num2str(inf_err_pos_S)])
        display(['max{|S(ij)-S*(ij)|: ij s.t. S*(ij)=0}: ',num2str(inf_err_0_S)])
    end
end