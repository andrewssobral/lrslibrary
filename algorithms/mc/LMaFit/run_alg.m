%{
load('dataset/trafficdb/traffic_patches.mat');
[M,m,n,p] = convert_video3d_to_2d(im2double(imgdb{100}));
out = run_algorithm('MC', 'LMaFit', M, [])
show_results(M.*out.Omega,out.L,out.S,out.O,p,m,n);
%}

MIdx = M(Idx);
rank = 10;
[X,Y] = lmafit_mc_adp(size(M,1),size(M,2),rank,Idx,MIdx,[]);
L = X*Y; % low-rank
S = (M - L); % sparse
