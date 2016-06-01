% Bayesian Robust PCA with Markov Dependency (Ding et al. 2011)
% process_video('RPCA', 'BRPCA-MD', 'dataset/demo.avi', 'output/demo_BRPCA-MD.avi');
%
K = 20;
Theta0 = InitialPara_random_MarkovDep(M,K);
[~,N] = size(M);
hyperpara.a0 = 1/K;
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
MCMCpara.nBurnin = 500;
MCMCpara.nCollect = 100;
%MCMCpara.nBurnin = 50;
%MCMCpara.nCollect = 10;
timerVal = tic;
out = Bayesian_RPCAmcmc_MarkovDep(M,Theta0,params.rows,params.cols,hyperpara,MCMCpara);
L = out.Lowrank_mean;
S = out.Sparse_mean;
