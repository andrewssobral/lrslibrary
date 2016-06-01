% ST | pROST | Robust PCA and subspace tracking from incomplete observations using L0-surrogates (Hage and Kleinsteuber, 2013)
% process_video('ST', 'pROST', 'dataset/demo.avi', 'output/demo_ST-pROST.avi');

L = robustpca_batch(M,2,'atansquare');
S = M - L;
