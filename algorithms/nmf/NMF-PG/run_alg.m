% NMF-PG: NMF solved by Projected Gradient
% process_video('NMF', 'NMF-PG', 'dataset/demo.avi', 'output/demo_NMF-PG.avi');
alg_path_aux = fullfile(lrs_conf.nmf_path,'NMF-DTU-Toolbox');
addpath(genpath(alg_path_aux));

M = sparse(M);

% cjlin: Alternative non-negative least squares using projected gradients.
[W, H] = nmf(M,2,'cjlin');

L = W * H;
S = M - L;

rmpath(genpath(alg_path_aux));
