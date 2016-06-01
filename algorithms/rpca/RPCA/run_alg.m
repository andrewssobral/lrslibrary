% RPCA (De la Torre and Black, 2001)
% process_video('RPCA', 'RPCA', 'dataset/demo.avi', 'output/demo_RPCA.avi');
sizeim = [params.rows params.cols];
[L,S] = RPCA(M,sizeim);
