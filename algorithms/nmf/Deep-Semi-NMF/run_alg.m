% Deep-Semi-NMF: Deep Semi Non-negative Matrix Factorization (Trigeorgis et al. 2014)
% process_video('NMF', 'Deep-Semi-NMF', 'dataset/demo.avi', 'output/demo_Deep-Semi-NMF.avi');
layers = 1;
[W,H] = deep_seminmf(M, layers);
W = cell2mat(W); H = cell2mat(H);
L = W * H;
S = M - L;
