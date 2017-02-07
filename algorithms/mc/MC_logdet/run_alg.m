%{
load('dataset/trafficdb/traffic_patches.mat');
[M,m,n,p] = convert_video3d_to_2d(im2double(imgdb{100}));
out = run_algorithm('MC', 'MC_logdet', M, [])
show_results(M.*out.Omega,out.L,out.S,out.O,p,m,n);
%}

mu = 0.1; % 0.1; 0.01;
kappa = 1.2; % 1.2; 1.0;
toler = 1e-10;
maxiter = 1000;
L = MC_LogDet_v3(M,Omega,mu,kappa,toler,maxiter);
S = M - L;
