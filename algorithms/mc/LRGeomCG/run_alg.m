%{
load('dataset/trafficdb/traffic_patches.mat');
[M,m,n,p] = convert_video3d_to_2d(im2double(imgdb{100}));
out = run_algorithm('MC', 'LRGeomCG', M, [])
show_results(M.*out.Omega,out.L,out.S,out.O,p,m,n);
%}

[L,S] = low_rank_matrix_completion(params);