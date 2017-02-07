%{
load('dataset/trafficdb/traffic_patches.mat');
[M,m,n,p] = convert_video3d_to_2d(im2double(imgdb{100}));
out = run_algorithm('MC', 'MC-NMF', M, [])
show_results(M.*out.Omega,out.L,out.S,out.O,p,m,n);
%}

MIdx = M(Idx);

[m,n] = size(M);
r = 2; t = 1;
esr = ceil(t*r); 

opts.tol = 1e-5;
opts.maxit = 500;
opts.print = 1;

[X,Y,Out] = mc_nmf(MIdx,Idx,esr,m,n,opts);
L = X*Y; % low-rank
S = (M - L); % sparse
