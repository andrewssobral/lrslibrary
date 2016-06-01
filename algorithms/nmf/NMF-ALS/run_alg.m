% NMF-ALS: NMF solved by Alternating Least Squares
% process_video('NMF', 'NMF-ALS', 'dataset/demo.avi', 'output/demo_NMF-ALS.avi');
alg_path_aux = fullfile(lrs_conf.nmf_path,'NMF-DTU-Toolbox');
addpath(genpath(alg_path_aux));

M = sparse(M);

% als: Alternating least squares.
[W, H] = nmf(M,1,'als');

L = W * H;
S = M - L;

rmpath(genpath(alg_path_aux));
