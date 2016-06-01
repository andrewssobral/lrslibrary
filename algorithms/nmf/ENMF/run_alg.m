% ENMF: Exact NMF (Gillis and Glineur, 2012)
% process_video('NMF', 'ENMF', 'dataset/demo.avi', 'output/demo_ENMF.avi');
rank = 1;
[H,W] = ExactNMF(M,rank,100);
L = (W' * H')';
S = M - L;
