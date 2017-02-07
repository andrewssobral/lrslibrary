%{
load('dataset/trafficdb/traffic_patches.mat');
[M,m,n,p] = convert_video3d_to_2d(im2double(imgdb{100}));
out = run_algorithm('MC', 'IALM-MC', M, [])
show_results(M.*out.Omega,out.L,out.S,out.O,p,m,n);
%}

MOmega = M.*Omega;
tol = 1e-4;
maxIter = 100;
A = inexact_alm_mc(MOmega,tol,maxIter);
L = A.U*A.V'; % low-rank
S = (M - L); % sparse
