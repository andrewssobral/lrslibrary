% PNMF: Probabilistic Non-negative Matrix Factorization
% process_video('NMF', 'PNMF', 'dataset/demo.avi', 'output/demo_PNMF.avi');
alg_path_aux = fullfile(lrs_conf.nmf_path,'NMF-DTU-Toolbox');
addpath(genpath(alg_path_aux));

M = sparse(M);

% prob: Probabilistic non-negative matrix factorization.
[W, H] = nmf(M,1,'prob');

L = W * H;
S = M - L;

rmpath(genpath(alg_path_aux));
