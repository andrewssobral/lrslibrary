%% LRSLibrary: A <Low-Rank and Sparse tools> Library for Background Modeling and Subtraction in Videos
close all; clear all; clc;

%% GUI
run_gui;

%% UTILS
video2mat('dataset/demo.avi', 'dataset/demo.mat');
video = load_video_file('dataset/demo.avi');
show_video(video);

%%% 2D demo
M = im2double(convert_video_to_2d(video));
show_2dvideo(M,video.height,video.width);

%%% 3D demo
video3d = im2double(convert_video_to_3d(video));
show_3dvideo(video3d);

%%% 3D tensor demo
add_tensor_libs;
T = convert_video_to_3dtensor(video);
slice3(double(T)), colormap('gray');
rem_tensor_libs;

%%% 4D demo
video4d = convert_video_to_4d(video);
video4d = crop_4dvideo(video4d,1,10);
video4d = resize_4dvideo(video4d,2);
show_4dvideo(video4d);

%% DEMO 01
load('dataset/trafficdb/traffic_patches.mat');
V = im2double(imgdb{100});
show_3dvideo(V);

%%% Matrix-based algorithms
[M, m, n, p] = convert_video3d_to_2d(V);
show_2dvideo(M,m,n);
% Robust PCA
results = process_matrix('RPCA', 'FPCP', M, []);
show_results(M,results.L,results.S,results.O,p,m,n);
% Low Rank Recovery
results = process_matrix('LRR', 'FastLADMAP', M, []);
show_results(M,results.L,results.S,results.O,p,m,n);
% Non-Negative Matrix Factorization
results = process_matrix('NMF', 'ManhNMF', M, []);
show_results(M,results.L,results.S,results.O,p,m,n);

%%% Tensor-based algorithms
add_tensor_libs;
T = tensor(V);
% Non-Negative Tensor Factorization
results = process_tensor('NTF', 'bcuNCP', T);
show_3dtensors(T,results.L,results.S,results.O);
% Tensor Decomposition
results = process_tensor('TD', 'Tucker-ALS', T);
show_3dtensors(T,results.L,results.S,results.O);
rem_tensor_libs;

%% DEMO 02

% Robust PCA
process_video('RPCA', 'FPCP', 'dataset/demo.avi', 'output/demo_FPCP.avi');
% Low Rank Recovery
process_video('LRR', 'FastLADMAP', 'dataset/demo.avi', 'output/demo_LRR-FastLADMAP.avi');
% Non-Negative Matrix Factorization
process_video('NMF', 'ManhNMF', 'dataset/demo.avi', 'output/demo_ManhNMF.avi');
% Non-Negative Tensor Factorization
process_video('NTF', 'bcuNCP', 'dataset/demo.avi', 'output/demo_bcuNCP.avi');
% Tensor Decomposition
process_video('TD', 'Tucker-ALS', 'dataset/demo.avi', 'output/demo_Tucker-ALS.avi');

%% DEMO 03 - For Large Videos (block by block)
close all; clear all; clc;
video = load_video_file('dataset/highway.avi');

%% show_video(video);
M_total = [];
L_total = [];
S_total = [];
O_total = [];
M = []; k = 1; k_max = 50;
nframes = 250; % video.nrFramesTotal;
for i = 1 : nframes
  %disp(['#frame ' num2str(i)]);
  frame = video.frames(i).cdata;
  if(size(frame,3) == 3)
    frame = rgb2gray(frame);
  end
  I = reshape(frame,[],1);
  M(:,k) = I;
  if(k == k_max || i == nframes)
    disp(['#last frame ' num2str(i)]);
    M = im2double(M);
    tic;
    results = process_matrix('RPCA', 'GoDec', M, []);
    %results = process_matrix('RPCA', 'IALM', M, []);
    %results = process_matrix('RPCA', 'FPCP', M, []);
    %results = process_matrix('LRR', 'FastLADMAP', M, []);
    %results = process_matrix('NMF', 'NMF-MU', M, []);
    toc
    M_total = [M_total M];
    L_total = [L_total results.L];
    S_total = [S_total results.S];
    O_total = [O_total results.O];
    displog('Displaying results...');
    show_results(M,results.L,results.S,results.O,size(M,2),video.height,video.width);
    M = []; k = 0;
    %break;
  end
  k = k + 1;
end
disp('Finished');

%%
show_results(M_total,L_total,S_total,O_total,size(M_total,2),video.height,video.width);

%%
convert_video2d_to_avi(S_total,size(S_total,2),video.height,video.width,'output/highway_S.avi');
convert_video2d_to_avi(L_total,size(L_total,2),video.height,video.width,'output/highway_L.avi');
convert_video2d_to_avi(O_total,size(O_total,2),video.height,video.width,'output/highway_O.avi');
