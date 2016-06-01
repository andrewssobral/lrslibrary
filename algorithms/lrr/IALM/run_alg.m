% LRR | IALM | Inexact ALM (Lin et al. 2009)
% process_video('LRR', 'IALM', 'dataset/demo.avi', 'output/demo_LRR-IALM.avi');

alg_path_aux = fullfile(lrs_conf.lrr_path,'ALM');
addpath(genpath(alg_path_aux));

A = mean(M,2);
lambda = 0.01;
[Z,E] = solve_lrr(M,A,lambda,0,1,1);
% M_hat = A*Z + E;
L = A*Z;
S = E;

rmpath(genpath(alg_path_aux));