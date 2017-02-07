%{
load('dataset/trafficdb/traffic_patches.mat');
[M,m,n,p] = convert_video3d_to_2d(im2double(imgdb{100}));
out = run_algorithm('MC', 'PG-RMC', M, [])
show_results(M.*out.Omega,out.L,out.S,out.O,p,m,n);
%}

avg = mean(mean(M, 1), 2);
M2 = M - avg;

[U_t, SV_t] = ncrmc(M2, Idx);
L = U_t * SV_t + avg; % low-rank
S = (M - L); % sparse
