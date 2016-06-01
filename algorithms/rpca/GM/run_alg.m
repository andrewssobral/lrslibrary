% GM (Hauberg et al. 2014)
% process_video('RPCA', 'GM', 'dataset/demo.avi', 'output/demo_GM.avi');
alg_path_aux = fullfile(lrs_conf.rpca_path,'GA');
addpath(genpath(alg_path_aux));

L = grassmann_median(M', 1);

L = nma_rescale(L,min(M(:)),max(M(:))); % 0.968 0.207
L = repmat(L,1,size(M,2));
S = M - L;

rmpath(genpath(alg_path_aux));
% show_2dvideo(M,m,n);
% show_2dvideo(L,m,n);
% show_2dvideo(S,m,n);
