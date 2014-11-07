function Output = Bayesian_RPCAmcmc_MarkovDep(X0,Theta0,Row,Col,hyperpara,MCMCpara)
%Bayesian_RPCAmcmc_withMarkovDep: Bayesian robust principle component analysis   
%implemented by MCMC, considering Markov dependency on the sparse term in time and space.
%
% Model: X0 = D*(Delta.*Z)*S + S2.*Z2 + E
% ----------------------------------------------------
%
%USAGE: Output = Bayesian_RPCAmcmc_withMarkovDep(X0,Theta0,Row,Col,hyperpara,MCMCpara)
%
%INPUT :  
%   X0: P x N, input data matrix. For the video application, every column is a frame of the video.
%  
%   Theta0: struct data including the initial parameters of the model. The parameters of the model are initialized as follows,
%         D = Theta0.D: P x K matrix as a dictionary for lowrank learning
%         S = Theta0.S: K x N coefficient matrix
%         Z = Theta0.Z: K x K diagnal binary matrix for rank learning
%         Delta = Theta0.Delta: K x K diagnal matrix for weighting Z.
%         Tao = Theta0.Tao: precision of Delta. Here we set Tao=1.
%         Pi = Theta0.Pi: K x 1 vector. Pi is the probability of Z=1.
%         
%         gamma_epsi = Theta0.gamma_epsi:  scalar, precision of the noise. Here we assumed the noise is stationary.
%         
%         S2 = Theta0.S2: P x N matrix for Sparse Component learning
%         Z2 = Theta0.Z2: P x N binary matrix for Sparse Component learning
%         gamma_s = Theta0.gamma_s: precision of the S2
%         Pi2 = Theta0.Pi2: P x N . Pi2 is the probability of Z2=1. 
%         Note that Pi2 here is a matrix, different from the model structure in function Bayesian_RPCAmcmc.
%
%   Row and Col: Row and column of the video frame.
%
%   hyperpara: MCMC hyperparameters
%     hyperpara.a0: scalar, hyperparameter 1 for Pi [1/K or 1/150 if K<150]
%       hyperpara.b0: scalar, hyperparameter 2 for Pi [1-hyperpara.a0]
%       hyperpara.c0: scalar, hyperparameter 1 for Tao precision [1e-6] 
%       hyperpara.d0: scalar, hyperparameter 2 for Tao precision [1e-6] 
%       hyperpara.e0: scalar, hyperparameter 1 for nosie precision [1e-6] 
%       hyperpara.f0: scalar, hyperparameter 2 for nosie precision [1e-6] 
%       hyperpara.g0: scalar, hyperparameter 1 for S2 precision [1e-6] 
%       hyperpara.h0: scalar, hyperparameter 2 for S2 precision [1e-6] 
%       hyperpara.alpha0: scalar,hyperparameter 1 for Pi2 [0.01*N]
%       hyperpara.beta0: scalar,hyperparameter 2 for Pi2 [0.99*N]
%       hyperpara.alpha1: scalar,hyperparameter 3 for Pi2 [0.99*N]
%       hyperpara.beta1: scalar,hyperparameter 4 for Pi2 [0.01*N]
%
%
%   MCMCpara: MCMC parameters
%       MCMCpara.nBurnin: scalar, number of burnin iterations [2000]
%       MCMCpara.nCollect: scalar, number of collected samples [1000]
%
%OUTPUT:
%   Output: struct data.
%         Output.Lowrank_mean: P x N matrix, mean of the Lowrank Component.
%         Output.Lowrank_std: P x N matrix, std of the Lowrank Component.
%         Output.Sparse_mean: P x N matrix, mean of the Sparse Component.
%         Output.Sparse_std: P x N matrix, std of the Sparse Component.
%         Output.Gamma_epsi_mean: scalar, mean of the noise precision.
%         Output.Gamma_epsi_std: scalar, std of the noise precision.
%         Output.rankL_mean: scalar, mean of the estimated rank of the lowrank Component.
%         Output.rankL_std: scalar, std of the estimated rank of the lowrank Component.
%         Output.NumSparse_mean: scalar, mean of the estimated number of the Sparse Component.
%         Output.NumSparse_std: scalar, std of the estimated number of the Sparse Component.


%--------------------------------------------------------------------------
% References:
% X. Ding, L. He and L. Carin, Bayesian Robust Principal Component Analysis, submitted to IEEE Trans. Image Processing (2010)  
%
% Xinghao Ding, Lihan He, ECE, Duke University
% Created: Apr. 12, 2010,
% Last change: Aug. 2, 2010. 
%--------------------------------------------------------------------------

% ---------------------
% check input arguments
% ---------------------

% P -- Dimension of the data vector
% N -- Number of samples
[P,N] = size(X0);
% K -- the largest possible rank
K = size(Theta0.D,2);

if nargin<6
    MCMCpara.nBurnin=100;
    MCMCpara.nCollect=100;
end
if nargin<5
    % Hyperparameters
    if K<150
        hyperpara.a0 = 1/150;
    else
        hyperpara.a0 = 1/K;
    end    
    hyperpara.b0 = 1-hyperpara.a0;
    hyperpara.c0 = 1e-6;
    hyperpara.d0 = 1e-6;
    hyperpara.e0 = 1e-6;
    hyperpara.f0 = 1e-6;
    hyperpara.g0 = 1e-6;
    hyperpara.h0 = 1e-6;    
    hyperpara.alpha0 = 0.01*N;
    hyperpara.beta0 = 0.99*N;
    hyperpara.alpha1 = 0.99*N;
    hyperpara.beta1 = 0.01*N; 
    
end

if isempty(hyperpara)
    % P -- Dimension of the data vector
    % N -- Number of samples
    [P,N] = size(X0);
    % K -- the largest possible rank
    K = size(Theta0.D,2);

    % Hyperparameters
    if K<150
        hyperpara.a0 = 1/150;
    else
        hyperpara.a0 = 1/K;
    end    
    hyperpara.b0 = 1-hyperpara.a0;
    hyperpara.c0 = 1e-6;
    hyperpara.d0 = 1e-6;
    hyperpara.e0 = 1e-6;
    hyperpara.f0 = 1e-6;
    hyperpara.g0 = 1e-6;
    hyperpara.h0 = 1e-6;
    hyperpara.alpha0 = 0.01*N;
    hyperpara.beta0 = 0.99*N;
    hyperpara.alpha1 = 0.99*N;
    hyperpara.beta1 = 0.01*N; 
end

if isempty(MCMCpara)
   MCMCpara.nBurnin=100;
    MCMCpara.nCollect=100;
end

a0 = hyperpara.a0;
b0 = hyperpara.b0;
c0 = hyperpara.c0;
d0 = hyperpara.d0;
e0 = hyperpara.e0;
f0 = hyperpara.f0;
g0 = hyperpara.g0;
h0 = hyperpara.h0;
alpha0 = hyperpara.alpha0;
beta0 =  hyperpara.beta0;
alpha1 = hyperpara.alpha1;
beta1 =  hyperpara.beta1; 


D = Theta0.D;
S = Theta0.S;
Z = Theta0.Z;
Delta = Theta0.Delta;
Tao = Theta0.Tao;
Pi = Theta0.Pi;
gamma_epsi = Theta0.gamma_epsi;
S2 = Theta0.S2;
Z2 = Theta0.Z2;
gamma_s = Theta0.gamma_s;
Pi2 = Theta0.Pi2;


Lowrank_Comp = zeros(P,N);
Sparse_Comp = zeros(P,N);


for iter=1:(MCMCpara.nBurnin+MCMCpara.nCollect)
    
    X = X0-D*diag(Delta.*Z)*S-S2.*Z2;  
    
    %---------------------------------------------------------------
    % Low-rank component sampling
    %---------------------------------------------------------------

    for k=1:K
        Tao(k) =1;  
        X = X+ Delta(k)*Z(k)*D(:,k)*S(k,:); 
        
        %Sample D----------------------------------------------------------
        sigma_Dk = 1/(gamma_epsi*(Delta(k)).^2*Z(k)*(S(k,:)*S(k,:)')+P);
        mu_Dk = gamma_epsi*sigma_Dk* (Delta(k))*Z(k)*X*S(k,:)';
        D(:,k) = mu_Dk + randn(P,1)*sqrt(sigma_Dk);
        %clear sigma_Dk mu_D
        %------------------------------------------------------------------   
        
        %Sample S----------------------------------------------------------
    
        Dk = D(:,k);
        sigS1 = 1/(1 + gamma_epsi*Z(k)*(Dk'*Dk)*(Delta(k)).^2);        
        muS1 = gamma_epsi*sigS1*Z(k)*Delta(k)*D(:,k)'*X;
        S(k,:) = randn(1,N)*sqrt(sigS1) + muS1;
        %------------------------------------------------------------------
        
        
        %Sample Delta(k)
        sig_Delta = 1/(Tao(k)+gamma_epsi*(Dk'*Dk)*(S(k,:)*S(k,:)'));
        mu_Delta = gamma_epsi*sig_Delta*D(:,k)'*X*S(k,:)';
        if Z(k) ==1
            %Delta(k) = randn(1)*sqrt(sig_Delta)+mu_Delta;
            Delta(k) = normrnd(mu_Delta,sqrt(sig_Delta));
        else
            %Delta(k) = randn(1);
            Delta(k) =  normrnd(0,sqrt(1/Tao(k)));
        end
        %------------------------------------------------------------------
       
%        %Sample Z---------------------------------------------------------
        Sk = S(k,:);
        Dk = D(:,k);
        temp =  - 0.5*gamma_epsi*(Delta(k)).^2*(Dk'*Dk)*(Sk*Sk') + gamma_epsi*Delta(k)*Dk'*(X*S(k,:)');          
        p1 = exp(temp)*Pi(k);      
        p0 = 1-Pi(k);        
        Z(k) = rand > p0/(p0+p1);
        %------------------------------------------------------------------
      
%       %Clappsed Gibbs Sample Z
%         tmpz = log(Pi(k)+eps) - log(1-Pi(k)+eps) + 0.5*log(sig_Delta) + (mu_Delta^2)/(2*sig_Delta) + 0.5*log(Tao(k));
%         if rand < 1/(1+exp(-tmpz))
%             Z(k) = 1; Delta(k) = normrnd(mu_Delta,sqrt(sig_Delta));
%         else
%             Z(k) = 0; Delta(k) =  normrnd(0,sqrt(1/Tao(k)));
%         end
        
        
        %sample Pi--------------------------------------------------
        ai = a0 + Z(k);
        bi = b0 + 1 - Z(k);
        Pi(k) =betarnd(ai,bi);
        %------------------------------------------------------------------
        
        %sample Tao
%         if Z(k) ==1
%             c1 =c0+1/2;
%             d1 =d0 + 0.5*(Delta(k)).^2;
%             Tao(k) = gamrnd(c1,1./d1);
%         else
%             Tao(k) = 1;
%         end
        %------------------------------------------------------------------
              
        X = X - Delta(k)*Z(k)*D(:,k)*S(k,:); 
    end  
    
    
    %---------------------------------------------------------------
    % Sparse component sampling
    %---------------------------------------------------------------

    X = X + S2.*Z2;  
    
    % Sample S2
    sig_S2 = 1./(gamma_s + gamma_epsi*Z2);
    mu_S2 = gamma_epsi*sig_S2.*Z2.*X;
    S2 = randn(P,N).*sqrt(sig_S2)+mu_S2;
    
    % Sample Z2    
    temp =  exp(-0.5*gamma_epsi*(S2.^2-2*S2.*X));    
%     p1 = repmat(Pi2,[1,N]).*temp;
%     p0 = repmat(1- Pi2,[1,N]);   
    p1 = Pi2.*temp;
    p0 = 1- Pi2;
    Z2 = rand(P,N)>p0./(p1+p0);
    
    X = X - S2.*Z2;  
    
    %sample gamma_s
    g1 = g0 + 0.5*P*N;
    tempS2 = S2.*S2;
    h1 = h0 + 0.5*sum(tempS2(:));
    gamma_s = gamrnd(g1,1/h1);
    
    
    %Sample Pi2
%     sumZ2 = sum(Z2,2);
%     a2 = a1 + sumZ2;
%     b2 = b1 + N -sumZ2;
%     Pi2 = betarnd(a2,b2);    
    
    %h = fspecial('average');
    h = [0.1,0.1,0.1;0.1,0.2,0.1;0.1,0.1,0.1];
    for n = 1:N
        temp1 = reshape(Z2(:,n),[Row,Col]);
        temp2 = imfilter(double(temp1),h)>0.6; 
        tempZ2(:,n) = reshape(temp2,[Row*Col,1]);
    end
    
    a2 = alpha0 +Z2(:,1);
    b2 = beta0 + 1 - Z2(:,1);
    Pi2(:,1) = betarnd(a2,b2);
    
    for t = 2:N
%         a2 = alpha0*(~Z2(:,t-1))+alpha1*(Z2(:,t-1)) +Z2(:,t);
%         b2 = beta0*(~Z2(:,t-1))+beta1*(Z2(:,t-1)) +1 - Z2(:,t);
%         a2 = alpha0*(~(Z2(:,t-1).*tempZ2(:,t-1)))+alpha1*(Z2(:,t-1).*tempZ2(:,t-1)) +Z2(:,t);
%         b2 = beta0*(~(Z2(:,t-1).*tempZ2(:,t-1)))+beta1*(Z2(:,t-1).*tempZ2(:,t-1)) +1 - Z2(:,t);
%         a2 = alpha0*(~(tempZ2(:,t-1)))+alpha1*(tempZ2(:,t-1)) +Z2(:,t);
%         b2 = beta0*(~(tempZ2(:,t-1)))+beta1*(tempZ2(:,t-1)) +1 - Z2(:,t);
          a2 = alpha0*(~(tempZ2(:,t)))+alpha1*(tempZ2(:,t)) +Z2(:,t);
          b2 = beta0*(~(tempZ2(:,t)))+beta1*(tempZ2(:,t)) +1 - Z2(:,t);
        
%         a2 = alpha0*(~(tempZ2(:,t-1).*tempZ2(:,t)))+alpha1*(tempZ2(:,t-1).*tempZ2(:,t)) +Z2(:,t);
%         b2 = beta0*(~(tempZ2(:,t-1).*tempZ2(:,t)))+beta1*(tempZ2(:,t-1).*tempZ2(:,t)) +1 - Z2(:,t);
        
%         a2 = alpha0*(~(Z2(:,t-1).*tempZ2(:,t)))+alpha1*(Z2(:,t-1).*tempZ2(:,t)) +Z2(:,t);
%         b2 = beta0*(~(Z2(:,t-1).*tempZ2(:,t)))+beta1*(Z2(:,t-1).*tempZ2(:,t)) +1 - Z2(:,t);
        
        Pi2(:,t) = betarnd(a2,b2);
    end
    
    
    %---------------------------------------------------------------
    % noise component sampling
    %---------------------------------------------------------------

    %Sample gamma_epsi
    e1 = e0 + 0.5*P*N;
    f1 = f0 + 0.5*sum(sum(X.^2));
    gamma_epsi = gamrnd(e1,1./f1);
    
    
    %------------------------
    
    % Collect samples
    if iter>MCMCpara.nBurnin
        ii = ceil(iter-MCMCpara.nBurnin);  
        
%         Lowrank_Comp(:,:,ii) = D*diag(Delta.*Z)*S;
%         Sparse_Comp(:,:,ii) = S2.*Z2;
        Lowrank_Comp = Lowrank_Comp + D*diag(Delta.*Z)*S;
        Sparse_Comp = Sparse_Comp + S2.*Z2;

        Gamma_epsi(ii) = gamma_epsi;
        tmpRank(ii) = length(find(Z~=0));
        %tmpNumSparse(ii) = length(find(Sparse_Comp(:,:,ii)~=0));
        %mse_rec(iter) = sum(sum((X0-Lowrank_Comp{ii}-Sparse_Comp{ii}).^2))/(P*N);
    end
    
end

% Output.Lowrank_mean = mean(Lowrank_Comp,3);
% Output.Lowrank_std = std(Lowrank_Comp,0,3);
% Output.Sparse_mean = mean(Sparse_Comp,3);
% Output.Sparse_std = std(Sparse_Comp,0,3);
Output.Lowrank_mean = Lowrank_Comp/MCMCpara.nCollect;
Output.Sparse_mean = Sparse_Comp/MCMCpara.nCollect;
Output.Gamma_epsi_mean = mean(Gamma_epsi);
Output.Gamma_epsi_std = std(Gamma_epsi);
Output.rankL_mean = mean(tmpRank);
Output.rankL_std = std(tmpRank);
% Output.NumSparse_mean = mean(tmpNumSparse);
% Output.NumSparse_std = std(tmpNumSparse);
end


