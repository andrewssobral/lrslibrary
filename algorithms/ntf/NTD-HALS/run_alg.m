% NTF | NTD-HALS | Non-negative Tucker Decomposition solved by Hierarchical ALS  (Zhou et al. 2012)
% process_video('NTF', 'NTD-HALS', 'dataset/demo.avi', 'output/demo_NTD-HALS.avi');

alg_path_aux = fullfile(lrs_conf.ntf_path,'lraNTD');
addpath(genpath(alg_path_aux));

R = [size(T,1) size(T,2) 2];
opts = struct('NumOfComp',R,'nlssolver','hals','maxiter',100,...
  'maxiniter',20,'tdalgFile','call_tucker_als_opts.mat');

T_hat = lraNTD_ANLS(T, opts);

L = double(tensor(T_hat));
S = double(T) - L;

rmpath(genpath(alg_path_aux));
