% NMF-MU: NMF solved by Multiplicative Updates
% process_video('NMF', 'NMF-MU', 'dataset/demo.avi', 'output/demo_NMF-MU.avi');
alg_path_aux = fullfile(lrs_conf.nmf_path,'NMF-DTU-Toolbox');
addpath(genpath(alg_path_aux));

M = sparse(M);

% mm: Multiplicative update method using euclidean distance measure.
[W, H] = nmf(M,1,'mm');

L = W * H;
S = M - L;

rmpath(genpath(alg_path_aux));
