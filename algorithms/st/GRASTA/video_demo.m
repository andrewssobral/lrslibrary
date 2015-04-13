%% This Demo intends to show that GRASTA can efficiently seperate the video into foreground / 
% background in realtime.
%
% Author: Jun He
% Email: hejun.zz@gmail.com
% 
% Papers used this demo code:
%[1] Jun He, Laura Balzano, and John C.S. Lui. Online robust subspace tracking from partial information. 
%    Preprint available at http://arxiv.org/pdf/1109.3827v2., 2011.
%
%[2] Jun He, Laura Balzano, and Arthur Szlam. Incremental gradient on the grassmannian for online foreground 
%    and background separation in subsampled video. In IEEE Conference on Computer Vision and Pattern Recognition 
%   (CVPR), June 2012
%

clear all; clc; close all;
grasta_path; % Add search path

%% Video parameters
% The following datasets can be download from 
% http://perception.i2r.a-star.edu.sg/bk_model/bk_index.html
%
DATASET         = 'lobby';

% Set your video path and groundtruth path like these
video_path    = '/Users/jhe/Documents/MATLAB/myCode/Data/CV';
gt_video_path = '/Users/jhe/Documents/MATLAB/myCode/Data/CV/GroundTruth';

if strcmp(DATASET,'lobby'),    
    VIDEOPATH = [video_path filesep 'lobby/'];
elseif strcmp(DATASET,'hall'),    
    VIDEOPATH = [video_path filesep 'hall/'];
elseif strcmp(DATASET, 'bootstrap'),    
    VIDEOPATH = [video_path filesep 'Bootstrap/'];
elseif strcmp(DATASET, 'escalator'),
    VIDEOPATH = [video_path filesep 'Escalator/'];
elseif strcmp(DATASET, 'campus'),    
    VIDEOPATH = [video_path filesep 'Campus/'];    
elseif strcmp(DATASET, 'curtain'),    
    VIDEOPATH = [video_path filesep 'Curtain/'];
elseif strcmp(DATASET, 'fountain'),    
    VIDEOPATH = [video_path filesep 'Fountain/'];
elseif strcmp(DATASET, 'watersurface'),    
    VIDEOPATH = [video_path filesep 'Watersurface/'];
elseif strcmp(DATASET, 'shopping'),    
    VIDEOPATH = [video_path filesep 'ShoppingMall/'];
end

% load the groundtruth frames
[frames_vec frame_names frame_names_gt] = groundTruthFrames(DATASET);
video_frame_name = videoFrames(DATASET);

%% GRASTA parameters
OPTIONS.RANK                = 5;  % the estimated low-rank
OPTIONS.rho                 = 1.8;    
OPTIONS.ITER_MAX            = 20; 
OPTIONS.ITER_MIN            = 20;    % the min iteration allowed for ADMM at the beginning

OPTIONS.USE_MEX             = 1;     % If you do not have the mex-version of Alg 2
                                     % please set Use_mex = 0.                                     

%% Initialize the rough subspace
OPTIONS.CONSTANT_STEP       = 0;   % use adaptive step-size to initialize the subspace
OPTIONS.MAX_LEVEL           = 20;
OPTIONS.MAX_MU              = 10000; % set max_mu large enough for initial subspace training
OPTIONS.MIN_MU              = 1;
FPS_ONLY                    = 1;    % 0:show the training video, 1: suppress the video demostration
TRAIN_FRAME                 = 1;    % 0¡Guse the first #training_size frames¡F 
                                    % 1: random select #training_size frames
                                    
max_cycles                  = 10;    % training cycles
training_size               = 100;   % random chose 50 frames as the training set
TRAINING_SAMPLING           = 0.3;   % Use how much information to train the first subspace.

t_start = tic;
[bgU, status, OPTS]  = bgtraining( VIDEOPATH, OPTIONS, max_cycles, TRAINING_SAMPLING, training_size,FPS_ONLY,TRAIN_FRAME);
toc(t_start);

%% Real-time background/foreground separation by Grasta
OPTS.MAX_ITER               = 20;
OPTIONS.CONSTANT_STEP       = 1e-2; % use the constant step-size
FPS_ONLY                    = 1;    % if you want to measure the FPS performance, please let FPS_ONLY=1
SAMPLING                    = 0.1;  % Use how much information to track the subspace.
thresh                      = 0.2;
MAX_FRAME                   = -1;   % -1 means seperating all the frames

OPTIONS.USE_MEX             = 1;
fprintf('Real-time seperation by Grasta with subsampling %.2f%%\n',100*SAMPLING);
fgs_grasta = bgfg_seperation_grasta( bgU, VIDEOPATH, SAMPLING ,status,OPTIONS, OPTS,thresh,FPS_ONLY, frame_names,MAX_FRAME);

%% Calculate the ROC curve

MAX_gt = 20; % use how many grounttruth frames to calculate the ROC curves [1--20] 

roc_grasta = calc_fgroc( gt_video_path, frames_vec, frame_names_gt, fgs_grasta,MAX_gt );

figure; plot(roc_grasta(:,1),roc_grasta(:,2),'r-o'); hold on;
xlabel('false positives'); ylabel('true positives');
legend('grasta'); title(['ROC curve of dataset -- '  DATASET]);
axis([0 0.2 0 .9]);grid on;

%% Make video -- grasta
OPTS.MAX_ITER               = 20;
OPTIONS.CONSTANT_STEP       = 1e-2; % use the constant step-size
FPS_ONLY                    = 1;    % if you want to measure the FPS performance, please let FPS_ONLY=1
SAMPLING                    = 0.2;  % Use how much information to track the subspace.
thresh                      = 0.2;
MAX_FRAME                   = -1;   % -1 means seperating all the frames
OPTIONS.USE_MEX             = 1;
fprintf('Seperating the whole video sequence by grasta...\n');
[video_grasta_fg,video_grasta_bg, vInfo] = bgfg_seperation_grasta( bgU, VIDEOPATH, SAMPLING ,status,OPTIONS, OPTS,thresh,FPS_ONLY, video_frame_name,MAX_FRAME);

%% Make video -- ReProCS(mod)
batch_size                  = 400; % training size for the initial clean background ReProCSmod by IALM
FPS_ONLY                    = 1;    % if you want to measure the FPS performance, please let FPS_ONLY=1
IALM_MAXITER                = 10;
MAX_FRAME                   = -1;
thresh                      = 0.2;
RESIZE                      = 1; % w = WIDTH/RESIZE, h=HEIGHT/RESIZE as ReProCSmod can not deal with large resolution video.
                                 % Cost too much memory! I have to resize
                                 % the video to low resolution...

fprintf('Seperating the whole video sequence by ReProCS(mod)...\n');
[video_reprocs_fg,video_reprocs_bg] = bgfg_seperation_ReProCSmod(VIDEOPATH,batch_size, IALM_MAXITER , thresh,FPS_ONLY , video_frame_name,MAX_FRAME,RESIZE);

%% Make video -- median filter
MEDIAN_BUF                  = 20;
FPS_ONLY                    = 1;
thresh                      = 0.2;
MAX_FRAME                   = -1;
fprintf('Seperating the whole video sequence by median filter...\n');
[video_mf_fg,video_mf_bg] = bgfg_seperation_mf(VIDEOPATH, MEDIAN_BUF, thresh,FPS_ONLY ,video_frame_name ,MAX_FRAME);