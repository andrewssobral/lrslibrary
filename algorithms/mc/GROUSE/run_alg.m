%{
load('dataset/trafficdb/traffic_patches.mat');
[M,m,n,p] = convert_video3d_to_2d(im2double(imgdb{100}));
out = run_algorithm('MC', 'GROUSE', M, [])
show_results(M.*out.Omega,out.L,out.S,out.O,p,m,n);
%}

maxrank = 1;
step_size = 0.1;
maxCycles = 100;
[Usg, Vsg, err_reg] = grouse(M,Omega,size(M,1),size(M,2),maxrank,step_size,maxCycles);
L = Usg*Vsg'; % low-rank
S = M - L; % sparse
