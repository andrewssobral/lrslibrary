% LRR | FastLADMAP | Fast LADMAP (Lin et al. 2011)
% process_video('LRR', 'FastLADMAP', 'dataset/demo.avi', 'output/demo_LRR-FastLADMAP.avi');

alg_path_aux = fullfile(lrs_conf.lrr_path,'ADM');
addpath(genpath(alg_path_aux));

lambda = 0.1;
rho = 1.9;
DEBUG = 1;
% X = XZ+E
[Z,E] = ladmp_lrr_fast(M,lambda,rho,DEBUG); clearvars -global A Xg eta M;
%M_hat = M*Z + E;
L = M*Z;
S = E; %S = M_hat - L;

rmpath(genpath(alg_path_aux));