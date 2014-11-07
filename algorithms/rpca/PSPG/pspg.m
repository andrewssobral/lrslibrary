function [X,S,out]=pspg(D,stdev,tol,optimal_X,optimal_S)

%************************************************************
% stopping criterion: |(X,S)-(Xp,Sp)|_F/(|(Xp,Sp)|_F+1)< tol*stdev
%************************************************************
% Written by Necdet Serhat Aybat, Penn State University
%  
% 1)first created on 28 August 2013.
%************************************************************

    OPTIONS.increaseK = 10;
    if nargin == 3
        output_flag = 0;
    elseif nargin == 5
        output_flag = 1;
    end
    %addpath(genpath('./'));
    nrm_D=norm(D,'fro');
    
    %tic
    [n1,n2] = size(D);
    n_max = max(n1,n2);
    n_min = min(n1,n2);
    
    deltabar = sqrt(n_min+sqrt(8*n_min))*stdev;
    
    norm_two = lansvd(D, 1, 'L');
    mu = 1.25/norm_two;     % this one can be tuned (1/smoothing parameter)
    mu_bar = mu * 1e8;
    coef = 1.5;            % this one can be tuned
    tau = 1;
    
    iter = 1; 
    flag=true;
    sv_array = [];
    rank_array = [];
    sv = round(n_min/99);
    
    X=0; Q=0;
    lambda = lambda_search(abs(D(:)), mu, deltabar, 1/sqrt(n_max));
    rho = (mu+lambda)/(sqrt(n_max)*mu*lambda);
        
    S = sign(D).*max(abs(D)-rho,0);
    X = lambda/(lambda+mu)*(D-S);
    Y = X;
    
    while(flag==true)
        
        OPTIONS.minSingValue = 1/mu;
        sv = round(0.99*sv);
        [U,Sigma,V] = lansvd(Y, sv, 'T', OPTIONS);
        Sigma = diag(Sigma);
        sv_array = [sv_array, length(Sigma)];
        sv = length(find(Sigma > 1/mu));
        rank_array = [rank_array, sv];

        Q = U(:, 1:sv) * diag(Sigma(1:sv) - 1/mu) * V(:, 1:sv)';
        
        %An alternative memory efficient way to compute Q
        %Q = bsxfun(@times,U(:, 1:sv),(diagSigma(1:sv) - 1/mu)') * V(:, 1:sv)';
        
        S_old = S;
        S = D - Q;
        lambda = lambda_search(abs(S(:)), mu, deltabar, 1/sqrt(n_max));
        rho = (mu+lambda)/(sqrt(n_max)*mu*lambda);
        
        S = sign(S).*max(abs(S)-rho,0);
        X_old = X;
        X = lambda/(lambda+mu)*(D-S)+ mu/(lambda+mu)*Q;
        
        stop_crit = sqrt(norm(X-X_old,'fro')^2+norm(S-S_old,'fro')^2)/(sqrt(norm(X_old,'fro')^2+norm(S_old,'fro')^2)+1);
        disp(['Iteration no: ',num2str(iter),' --> |(X,S)-(Xp,Sp)|_F/(1+|(Xp,Sp)|_F)=', num2str(stop_crit)])
        converged = stop_crit <= tol*stdev & iter>10;
        
        if converged
            %total_time = toc
            disp([num2str(iter),' itererations done'])
            flag=false;
        end
        tau_p = tau;
        tau = (1+sqrt(1+4*tau^2))/2;
        c = (tau_p-1)/tau; 
        Y = (c+1)*X-c*X_old;
        mu = min(mu*coef, mu_bar);
        iter=iter+1;
    end
    
    disp([num2str(iter-1),' total partial svd number'])
    if output_flag == 0
        [X,S] = treshold(X,D,stdev);
    end
    noise_treshold =1e-5;
    Sigma_X = svd(X,'econ');
    o1 = sum(Sigma_X);
    o2 = norm(S(:),1);
    o3 = norm(X+S-D,'fro');

    if output_flag == 0
        
        out.rel_err_X = norm(X-D,'fro')/nrm_D;
        out.rel_err_S = norm(S-D,'fro')/nrm_D;
        out.rel_err_noise = o3/nrm_D;
        out.rank = sum(Sigma_X>3*stdev);
        
        display(' ')
        display('**********************************************************')
        display(['Nuclear Norm: ', num2str(o1), '  L1 Norm: ', num2str(o2)])
        display(['|X+S-D|_F: ', num2str(o3)])
        display(['|D-X|/|D|: ',num2str(out.rel_err_X)])
        display(['|D-S|/|D|: ',num2str(out.rel_err_S)])
        display(['|X+S-D|/|D|: ',num2str(out.rel_err_noise)])
        
    else       
        Sigma_opt = svd(optimal_X, 'econ');
        nuc_opt = sum(Sigma_opt);
        l1_opt = norm(optimal_S(:),1);
        out.rel_err_X = norm(X-optimal_X,'fro')/norm(optimal_X,'fro');
        out.rel_err_S = norm(S-optimal_S,'fro')/norm(optimal_S,'fro');
        nuc_norm_rel_err = abs(nuc_opt-o1)/nuc_opt;
        L1_rel_err = abs(l1_opt-o2)/l1_opt;
        out.rel_err_noise = o3/nrm_D;
        out.rank = sum(Sigma_X>3*stdev);

        display(' ')
        display('**********************************************************')
        display(['Nuclear Norm(X*): ', num2str(nuc_opt),'  Nuclear Norm: ', num2str(o1)])
        display(['L1 Norm(S*): ', num2str(l1_opt),'  L1 Norm: ', num2str(o2)])
        display(['  |X+S-D|_F: ', num2str(o3)])
        display(['Relative Error X: ',num2str(out.rel_err_X)])
        display(['Relative Error s: ',num2str(out.rel_err_S)])
        
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
        
        display(' ')
        display(['|X+S-D|/|D|: ',num2str(out.rel_err_noise)])
    end   
end