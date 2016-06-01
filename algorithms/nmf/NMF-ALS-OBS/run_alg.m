% NMF-ALS-OBS: NMF solved by Alternating Least Squares with Optimal Brain Surgeon
% process_video('NMF', 'NMF-ALS-OBS', 'dataset/demo.avi', 'output/demo_NMF-ALS-OBS.avi');
alg_path_aux = fullfile(lrs_conf.nmf_path,'NMF-DTU-Toolbox');
addpath(genpath(alg_path_aux));

M = sparse(M);

% alsobs: Alternating least squares with optimal brain surgeon.
[W, H] = nmf(M,1,'alsobs');

L = W * H;
S = M - L;

rmpath(genpath(alg_path_aux));
