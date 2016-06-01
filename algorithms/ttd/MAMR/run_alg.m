% TTD | MAMR | Motion-Assisted Matrix Restoration (Ye et al. 2015)
% process_video('TTD', 'MAMR', 'dataset/demo.avi', 'output/demo_MAMR.avi');

%{
clc;
lrs_load_conf;
load(fullfile(lrs_conf.lrs_dir,'dataset','trafficdb','traffic_patches.mat'));
V = im2double(imgdb{100});
[M,m,n,p] = convert_video3d_to_2d(V);
%}

%L = zeros(size(M));
%S = zeros(size(M));
E = zeros(size(M));

%CM = ones(size(M));
%CM = zeros(size(M));
[numr,numc] = size(M);
CM = randi([0 1],numr,numc); % simulated confidence map (binary matrix)

Omega = find(CM ~= 0);
[I, J] = ind2sub([numr numc],Omega);

% Motion-Assisted Matrix Restoration (Ye et al. 2015)
lambda = 1;
[L,S,iter1] = core_MAMR(M,lambda,I,J);

M_hat = L + S + E;
error = norm(M_hat(:)-M(:))/norm(M(:));
disp(['Error: ' num2str(error)]);

%{
show_2dvideo(M,m,n);
show_2dvideo(M.*CM,m,n);
show_2dvideo(L,m,n);
show_2dvideo(S,m,n);
show_2dvideo(Z,m,n);
%}