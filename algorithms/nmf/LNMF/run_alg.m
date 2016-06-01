% LNMF: Spatially Localized NMF (Li et al. 2001)
% process_video('NMF', 'LNMF', 'dataset/demo.avi', 'output/demo_LNMF.avi');
rank = 1;
option.verbose = 1;
[W,H] = LNMF(M,rank,option);
L = W * H;
S = M - L;
