% Semi-NMF: Semi Non-negative Matrix Factorization (Trigeorgis et al. 2014)
% process_video('NMF', 'Semi-NMF', 'dataset/demo.avi', 'output/demo_Semi-NMF.avi');
rank = 10;
[W,H] = seminmf(M, rank);
L = W * H;
S = M - L;
