%{
load('dataset/trafficdb/traffic_patches.mat');
[M,m,n,p] = convert_video3d_to_2d(im2double(imgdb{100}));
out = run_algorithm('MC', 'OptSpace', M, [])
show_results(M.*out.Omega,out.L,out.S,out.O,p,m,n);
%}

MOmega = M.*Omega;
tol = 1e-8;
[X,S,Y,dist] = OptSpace(sparse(MOmega),[],20,tol);
L = X*S*Y'; % low-rank
S = M - L; % sparse
