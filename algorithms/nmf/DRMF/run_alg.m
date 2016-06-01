% DRMF: Direct Robust Matrix Factorization (Xiong et al. 2011)
% process_video('NMF', 'DRMF', 'dataset/demo.avi', 'output/demo_DRMF.avi');

%%% initialization
lambda = 1/sqrt(max(size(M)));
L_rpca = inexact_alm_rpca(M, lambda, 1e-5, 10);
sv = svdex(L_rpca);
rk = EffRank(sv, 0.999);
%%% run
options.init = L_rpca;
[L,S] = DRMF(M, rk, 0.1, options);
S = full(S);
%%% end
