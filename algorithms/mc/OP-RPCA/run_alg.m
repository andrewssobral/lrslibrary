% MC | OP-RPCA | Robust PCA via Outlier Pursuit (Xu et al. 2012)
% process_video('MC', 'OP-RPCA', 'dataset/demo.avi', 'output/demo_OP-RPCA.avi');

%{
clear, clc;
load('dataset/trafficdb/traffic_patches.mat');
V = im2double(imgdb{100});
[M,m,n,p] = convert_video3d_to_2d(V);
%}

%[numr,numc] = size(M);
%Idx = randi([0 1],numr,numc);
Idx = ones(size(M));

lambda = 0.35; % [0,1] 0.5(+lowrank) 0.4 0.3 0.2(+sparse)
[L,S] = mr_pca_part(M,Idx,lambda);

%{
show_2dvideo(M,m,n);
show_2dvideo(L,m,n);
show_2dvideo(S,m,n);
show_2dvideo(hard_threshold(S),m,n);
%}