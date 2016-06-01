% ST | GRASTA | Grassmannian Robust Adaptive Subspace Tracking Algorithm (He et al. 2012)
% process_video('ST', 'GRASTA', 'dataset/demo.avi', 'output/demo_ST-GRASTA.avi');

% clear, clc;
% load('dataset/trafficdb/traffic_patches.mat');
% V = im2double(imgdb{100});
% [M,m,n,p] = convert_video3d_to_2d(V);
% I = reshape(M(:,1),m,n); % imshow(I);

%%% predefine your video parameters
%VIDEO_ROWS  =  320; % Row of each frame
%VIDEO_COLS  =  240; % Column of each fram
%DIM         = VIDEO_ROWS * VIDEO_COLS;
DIM = size(M,1);

%%% GRASTA parameters
subsampling                 = 1;   % how much patial information will be used in your application
OPTIONS.RANK                = 1;     % the estimated rank
OPTIONS.rho                 = 1.8;
OPTIONS.MAX_MU              = 10000; % set max_mu large enough for initial subspace training
OPTIONS.MIN_MU              = 1;
OPTIONS.ITER_MAX            = 20;
OPTIONS.DIM_M               = DIM;   % your data's dimension
OPTIONS.USE_MEX             = 0;     % If you do not have the mex-version of Alg 2
% please set Use_mex = 0.
OPTS                        = struct(); % initiate a empty struct for OPTS
status.init                 = 0;        % status of grasta at each iteration
U_hat                       = zeros(1); % initiate a zero U_hat at beginning

%%% Initial subspace training [optional]
% You may use some frames to train the initial subspace.
OPTIONS.CONSTANT_STEP       = 0; % use adaptive step-size for initial subspace training
max_cycles                  = 30;
training_frames             = 10;
for outiter = 1:max_cycles
  frame_order = randperm(training_frames);
  for i = 1:training_frames
    % prepare the training frame
    %I = imread(fname);
    %I = double(rgb2gray(I));
    %I = I/max(max(I));
    I = M(:,frame_order(i));
    
    % random subsampling the frame I
    Z = round(subsampling * DIM);
    rp = randperm(DIM);
    idx = rp(1:Z)';
    I_Omega = I(idx);
    
    % training the background
    [U_hat,status,OPTS] = grasta_stream(I_Omega,idx,U_hat,status,OPTIONS,OPTS);
    %fprintf('Training %d/%d ...\n',outiter,max_cycles);
  end
end

%%% Real-time background/foreground separation
OPTIONS.CONSTANT_STEP = 1e-2; % use small constant step-size
nframes = size(M,2);
for i = 1:nframes
  % prepare the image I whether it is saved as file or caputured by
  % camera
  % I = imread(fname);
  %I = double(rgb2gray(I));
  %I = I/max(max(I));
  I = M(:,i);
  
  % random subsampling the frame I
  Z = round(subsampling * DIM);
  rp = randperm(DIM);
  idx = rp(1:Z)';
  I_Omega = I(idx);
  
  % tracking the background
  [U_hat,status,OPTS] = grasta_stream(I_Omega,idx,U_hat,status,OPTIONS,OPTS);
  
  % bg_img is the background
  L_hat = U_hat * status.w * status.SCALE;
  %bg_img = reshape(U_hat * status.w * status.SCALE, VIDEO_ROWS,VIDEO_COLS);
  
  % s_img is the separated foreground
  %s_hat = I(:) - U_hat * status.w * status.SCALE;
  %s_img = reshape(s_hat,VIDEO_ROWS,VIDEO_COLS);
  S_hat = I - L_hat;
  
  L(:,i) = L_hat;
  S(:,i) = S_hat;
  
  %fprintf('Processing %d/%d ...\n',i,nframes);
end

%%
% show_results(M,L,S,hard_threshold(S),p,m,n);