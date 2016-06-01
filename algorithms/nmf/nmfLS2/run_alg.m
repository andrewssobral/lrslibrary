% nmfLS2: Non-negative Matrix Factorization with sparse matrix (Ji and Eisenstein, 2013)
% process_video('NMF', 'nmfLS2', 'dataset/demo.avi', 'output/demo_nmfLS2.avi');
rank = 1;
[W,H] = nmfLS2(M, rank);
L = W * H;
S = M - L;
