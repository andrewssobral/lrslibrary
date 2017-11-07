function [BG, FG, L_hat, S_hat, T_hat, t_hat, ...
    P_track_full, P_track_new, T_calc]= ...
    ReProCS(M1, P_init, mu, ev_thresh, alpha, K)
%%%This is the ReProCS algorithm under the latest subspace change and 
%%%support model. This is the main function for the real video Background-
%%%Foreground separation problem. 

%This folder contains the code accompanying pre-print.
%
%[1] "New Results for Provable Dynamic Robust PCA", Praneeth Narayanamurthy and Namrata Vaswani, arXiv:1705.08948, 2017.
%
%If you use this code please also cite the following papers
%[2] "An online algorithm for separating sparse  and low-dimensional signal sequences from their sum", Han Guo, Chenlu Qiu, and Namrata Vaswani, IEEE Trans. Sig. Proc., 2014.
%[3] "Recursive Robust PCA or Recursive Sparse Recovery in Large but Structure Noise", Chenlu Qiu, Namrata Vaswani, Brain Lois, and Leslie Hogben, IEEE Trans. Info. Theory., 2014.
%[4] "Real-time Robust Principal Components' Pursuit", Chenlu Qiu, and Namrata Vaswani, Allerton, 2010.


%%%                          Inputs                         %%%
%%%     M - measurement matrix                              %%%
%%%     ev_thres - threshold for subspace change detection  %%%
%%%     P_init - an initial estimate of the subspace        %%%
%%%     t_train - the dimension of the training data        %%%


%%%                       Algorithm parameters              %%%
%%%     alpha - frame length                                %%%
%%%     mu - column-averages of data                        %%%
%%%     K - number of projection PCA steps                  %%%
%%%     omega - threshold for non-zero value in S           %%%


%%%                          Outputs                        %%%
%%%     BG - Estimated Background                           %%%
%%%     FG - Estimated Foreground                           %%%
%%%     L_hat - estimate of the low rank matrix             %%%
%%%     P_hat - estimate of the subspace in which L lies    %%%
%%%     S_hat - estimate of the sparse signal               %%%
%%%     t_hat - estimate of subspace change times           %%%

%% Initializations

P_hat_old = P_init;
P_hat_new = [];
P_hat = [P_hat_old, P_hat_new];

[n, t_max] = size(M1);
T_hat = zeros(n, t_max);
S_hat = zeros(n, t_max);
L_hat = zeros(n, t_max);

FG = zeros(n, t_max);
BG = zeros(n, t_max);

t_hat = [];
M = M1 - repmat(mu, 1, t_max);
k = 0;
cnt = 1;
ph = 0;     %ph - 0 => detect, 1 => ppca
opts.delta = 0.4;

phi_t = speye(n) - P_hat * P_hat';
%% Main Algorithm Loop
for ii = 2 : t_max
        %% Estimate support
        Atf.times = @(x) x - (P_hat * (P_hat' * x));
        Atf.trans = @(y) y - (P_hat * (P_hat' * y));
        phi.times = @(x) x - (P_hat_old * (P_hat_old' * x));
        y_t = Atf.times(M(:, ii));
        opts.tol   = 1e-3; 
        opts.print = 0;
        
        opts.delta = norm(Atf.times(L_hat(:, ii - 1)));
        
        x_t_hat_cs = yall1(Atf, y_t, opts); 
        omega = sqrt(M(:, ii)' * M(:, ii) / n);
        
        t_hat_temp = find(abs(x_t_hat_cs) > omega);
        T_hat(t_hat_temp, ii) = 255;
        
        LS.times = @(x) phi(:, t_hat_temp) * x;
        LS.trans = @(y) phi(:, t_hat_temp)' * x;
        
        %% Estimate signal components
        
        S_hat(t_hat_temp, ii) = cgls(phi_t(:, t_hat_temp), y_t, 0, 1e-3);
        L_hat(:, ii) = M(:, ii) - S_hat(:, ii);
        FG(t_hat_temp, ii) = M1(t_hat_temp, ii);
        BG(:, ii) = L_hat(:, ii) + mu;
        %% Subspace update
        if(~mod(ii - 1 , alpha))
            MM = (1 / sqrt(alpha)) * phi.times(L_hat(:, ii - alpha + 1 : ii));
            
            if(~ph)     %%detect phase
                aa = svds(MM, 1);
                if(aa >= sqrt(ev_thresh))
                    ph = 1;
                    t_hat = [t_hat, ii];
                    k = 0;
                end
            else        %%ppca phase
                P_hat_new = proj_PCA_thresh(MM, sqrt(alpha *ev_thresh));
                P_hat = [P_hat_old, P_hat_new];
                phi_t = speye(n) - P_hat * P_hat';
                k = k + 1;
                
                if(k==K + 1)
                    P_hat = [P_hat_old, P_hat_new];
                    P_hat_old = P_hat;
                    P_hat_new = [];
                    k = 1;
                    ph = 0;
                    phi_t = speye(n) - P_hat * P_hat';
                end
            end
        end
    
    %% Return subspace estimates
    if((ii == 0) || (~(mod(ii - 1, alpha))))
        P_track_full{cnt} = P_hat;
        P_track_new{cnt} = P_hat_new;
        T_calc(cnt) = ii;
        cnt = cnt + 1;
    end
end
%L_hat = L_hat + repmat(mu, n, 1);
end


