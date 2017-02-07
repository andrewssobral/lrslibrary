%{
load('dataset/trafficdb/traffic_patches.mat');
[M,m,n,p] = convert_video3d_to_2d(im2double(imgdb{100}));
out = run_algorithm('MC', 'SVP', M, [])
show_results(M.*out.Omega,out.L,out.S,out.O,p,m,n);
%}

warning('off','all');
MIdx = M(Idx);
rank = 2;
[U,S,V] = svp(Idx,MIdx,size(M,1),size(M,2),rank);
L = U*diag(S)*V'; % low-rank
S = (M - L); % sparse
warning('on','all');
