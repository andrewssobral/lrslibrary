% TFOCS (Becker et al. 2011)
% process_video('RPCA', 'TFOCS-IC', 'dataset/demo.avi', 'output/demo_TFOCS-IC.avi');

alg_path_aux = fullfile(lrs_conf.rpca_path,'TFOCS');
addpath(genpath(alg_path_aux));

L = tfocs_interface(M, 2);
S = M - L;

rmpath(genpath(alg_path_aux));
