% TFOCS (Becker et al. 2011)
% process_video('RPCA', 'TFOCS-EC', 'dataset/demo.avi', 'output/demo_TFOCS-EC.avi');

alg_path_aux = fullfile(lrs_conf.rpca_path,'TFOCS');
addpath(genpath(alg_path_aux));

L = tfocs_interface(M, 1);
S = M - L;

rmpath(genpath(alg_path_aux));
