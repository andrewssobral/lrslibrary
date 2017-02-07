%{
load('dataset/trafficdb/traffic_patches.mat');
[M,m,n,p] = convert_video3d_to_2d(im2double(imgdb{100}));
out = run_algorithm('MC', 'FPC', M, [])
show_results(M.*out.Omega,out.L,out.S,out.O,p,m,n);
%}

maxiter = 500;
mu_final = .01; tol = 1e-6;
MIdx = M(Idx);

fprintf('\nSolving by FPC...\n');
[U,S,V,numiter] = FPC(size(M),Idx,MIdx,mu_final,maxiter,tol);

L = U*S*V'; % low-rank
S = M - L; % sparse
