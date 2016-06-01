% LRR | LADMAP | Linearized ADM with Adaptive Penalty (Lin et al. 2011)
% process_video('LRR', 'LADMAP', 'dataset/demo.avi', 'output/demo_LRR-LADMAP.avi');

alg_path_aux = fullfile(lrs_conf.lrr_path,'ADM');
addpath(genpath(alg_path_aux));

lambda = 0.1;
rho = 1.9;
DEBUG = 1;
% X = XZ+E
[Z,E] = ladmp_lrr(M,lambda,rho,DEBUG); clearvars -global M;
%M_hat = M*Z + E;
L = M*Z;
S = E; %S = M_hat - L;

rmpath(genpath(alg_path_aux));