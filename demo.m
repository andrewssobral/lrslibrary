%% LRSLibrary: A library of low-rank and sparse tools for background/foreground separation in videos
close, clear, clc;
% restoredefaultpath;

%% First run the setup script
lrs_setup; % or run('C:/GitHub/lrslibrary/lrs_setup')

%% LRS GUI (graphical user interface)
lrs_gui;

%% Load configuration (for demos)
lrs_load_conf;

%% Auxiliary functions

%%% Load video
input_avi = fullfile(lrs_conf.lrs_dir,'dataset','demo.avi');
video = load_video_file(input_avi);
show_video(video);

%%% Convert video to MAT (MATLAB format)
output_mat = fullfile(lrs_conf.lrs_dir,'dataset','demo.mat');
video2mat(input_avi, output_mat);

%%% Convert video to matrix representation (2D array)
M = im2double(convert_video_to_2d(video));
show_2dvideo(M,video.height,video.width);

%%% Convert video to 3D array
V = im2double(convert_video_to_3d(video));
show_3dvideo(V);

%%% Convert video to 3D tensor
T = convert_video_to_3dtensor(video);
tensorlab.slice3(double(T)), colormap('gray'); % Visualize a third-order tensor with slices.

%%% Convert video to 4D array/tensor
video4d = convert_video_to_4d(video);
video4d = crop_4dvideo(video4d,1,10);
video4d = resize_4dvideo(video4d,2);
show_4dvideo(video4d);

%% DEMO 01 (process matrix/tensor data)
load(fullfile(lrs_conf.lrs_dir,'dataset','trafficdb','traffic_patches.mat'));
V = im2double(imgdb{100});
show_3dvideo(V);

%%% For matrix-based algorithms
[M,m,n,p] = convert_video3d_to_2d(V);
show_2dvideo(M,m,n);

% Robust PCA
out = run_algorithm('RPCA', 'FPCP', M);
% Subspace Tracking
out = run_algorithm('ST', 'GRASTA', M);
% Matrix Completion (sample randomly 50% of the observations)
out = run_algorithm('MC', 'GROUSE', M);
% Low Rank Recovery
out = run_algorithm('LRR', 'FastLADMAP', M);
% Three-Term Decomposition
out = run_algorithm('TTD', '3WD', M);
% Non-Negative Matrix Factorization
out = run_algorithm('NMF', 'ManhNMF', M);

% Show results
show_results(M,out.L,out.S,out.O,p,m,n);

%%% Tensor-based algorithms
T = tensor(V);

% Non-Negative Tensor Factorization
out = run_algorithm('NTF', 'bcuNCP', T);
% Tensor Decomposition
out = run_algorithm('TD', 'Tucker-ALS', T);

% Show results
show_3dtensors(T,out.L,out.S,out.O);

%% DEMO 02 (data completion)
load(fullfile(lrs_conf.lrs_dir,'dataset','trafficdb','traffic_patches.mat'));
V = im2double(imgdb{100});
show_3dvideo(V);

%%% For matrix-based algorithms
[M,m,n,p] = convert_video3d_to_2d(V);
show_2dvideo(M,m,n);

%%% Define observed entries
obs = 0.5; % Percentual [0...1]
[params.Idx, params.Omega] = subsampling(M, obs);

%%% Matrix completion
out = run_algorithm('MC', 'GROUSE', M, params);

%%% Show results
show_results(M.*out.Omega,out.L,out.S,out.O,p,m,n);

%% DEMO 03 (process video data)

% Load video
input_avi = fullfile(lrs_conf.lrs_dir,'dataset','demo.avi');
output_avi = fullfile(lrs_conf.lrs_dir,'output','output.avi');

% Robust PCA
process_video('RPCA', 'FPCP', input_avi, output_avi);
% Subspace Tracking
process_video('ST', 'GRASTA', input_avi, output_avi);
% Matrix Completion
process_video('MC', 'GROUSE', input_avi, output_avi);
% Low Rank Recovery
process_video('LRR', 'FastLADMAP', input_avi, output_avi);
% Three-Term Decomposition
process_video('TTD', '3WD', input_avi, output_avi);
% Non-Negative Matrix Factorization
process_video('NMF', 'ManhNMF', input_avi, output_avi);
% Non-Negative Tensor Factorization
process_video('NTF', 'bcuNCP', input_avi, output_avi);
% Tensor Decomposition
process_video('TD', 'Tucker-ALS', input_avi, output_avi);

%% DEMO 04 - Process large videos (block by block)
input_avi = fullfile(lrs_conf.lrs_dir,'dataset','highway.mpg');
video = load_video_file(input_avi);
% show_video(video);

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
    out = run_algorithm('RPCA', 'GoDec', M, []);
    %results = run_algorithm('RPCA', 'IALM', M, []);
    %results = run_algorithm('RPCA', 'FPCP', M, []);
    %results = run_algorithm('LRR', 'FastLADMAP', M, []);
    %results = run_algorithm('NMF', 'NMF-MU', M, []);
    toc
    M_total = [M_total M];
    L_total = [L_total out.L];
    S_total = [S_total out.S];
    O_total = [O_total out.O];
    displog('Displaying results...');
    show_results(M,out.L,out.S,out.O,size(M,2),video.height,video.width);
    M = []; k = 0;
    %break;
  end
  k = k + 1;
end
disp('Finished');

%% Show results
show_results(M_total,L_total,S_total,O_total,size(M_total,2),video.height,video.width);

%% Convert 2D matrix to AVI
convert_video2d_to_avi(S_total,size(S_total,2),video.height,video.width,'output/highway_S.avi');
convert_video2d_to_avi(L_total,size(L_total,2),video.height,video.width,'output/highway_L.avi');
convert_video2d_to_avi(O_total,size(O_total,2),video.height,video.width,'output/highway_O.avi');
