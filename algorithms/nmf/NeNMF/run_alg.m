% NeNMF: NMF via Nesterov's Optimal Gradient Method (Guan et al. 2012)
% process_video('NMF', 'NeNMF', 'dataset/demo.avi', 'output/demo_NeNMF.avi');
rank = 1;
[W,H] = NeNMF(M,rank);
L = W * H;
S = M - L;
