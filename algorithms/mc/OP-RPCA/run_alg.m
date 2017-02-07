%{
load('dataset/trafficdb/traffic_patches.mat');
[M,m,n,p] = convert_video3d_to_2d(im2double(imgdb{100}));
out = run_algorithm('MC', 'OP-RPCA', M, [])
show_results(M.*out.Omega,out.L,out.S,out.O,p,m,n);
%}

lambda = 0.35; % [0,1] 0.5(+lowrank) 0.4 0.3 0.2(+sparse)
[L,S] = mr_pca_part(M,Omega,lambda);
