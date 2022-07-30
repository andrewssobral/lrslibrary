%%%Wrapper for foreground-background separation 
%This folder contains the code accompanying pre-print.
%
%[1] "Provable Dynamic Robust PCA or Robust Subspace Tracking", Praneeth Narayanamurthy and Namrata Vaswani, arXiv:1705.08948, 2017.
%
%If you use this code please also cite the following papers
%[2] "An online algorithm for separating sparse  and low-dimensional signal sequences from their sum", Han Guo, Chenlu Qiu, and Namrata Vaswani, IEEE Trans. Sig. Proc., 2014.
%[3] "Recursive Robust PCA or Recursive Sparse Recovery in Large but Structure Noise", Chenlu Qiu, Namrata Vaswani, Brain Lois, and Leslie Hogben, IEEE Trans. Info. Theory., 2014.
%[4] "Real-time Robust Principal Components' Pursuit", Chenlu Qiu, and Namrata Vaswani, Allerton, 2010.


%%Read video
%clear;
%clc;
%close all

%addpath('YALL1_v1.4/');
%addpath('PROPACK/');

%load('Data/Curtain.mat');
%I = M;


%% Training data processing
%%option 1 -- init using batch RPCA
t_train = 10; %400;
TrainData = M(:, 1 : t_train);
rank_init = 5;

L_hat_init = ncrpca(TrainData, rank_init);

mu = mean(L_hat_init, 2);


[Utemp, Stemp, ~] = svd(1 / sqrt(t_train) * (L_hat_init - ...
    repmat(mu, 1, t_train)));
ss1 = diag(Stemp);
L_init = Utemp(:, 1 : rank_init);

%% option 2 -- init using outlier free data
% mu = mean(DataTrain, 2);
% t_train = size(DataTrain, 2);
% [Utemp, Stemp, ~] = svd(1 / sqrt(t_train) * ...
%     (DataTrain - repmat(mu, 1, t_train)));
% ss1 = diag(Stemp);
% b = 0.95;
% rank_init = min(find(cumsum(ss1.^2) >= b * sum(ss1.^2)));
% L_init = Utemp(:, 1 : rank_init);


fprintf('Initialized\n');
theta_thresh = 20 * pi / 180;
ev_thresh = 0.1 * ss1(rank_init) * sin(theta_thresh)^2;


%%Main online robust PCA algorithm section
%% Call to online RPCA function

K = 3;
alpha = 60;


tic
fprintf('alpha = %d\tK = %d\n', alpha, K);
%[BG, FG, L_hat, S_hat, T_hat, t_hat, P_track_full, P_track_new] ...
%    = ReProCS(M(:, t_train + 1 : end), ...
%    L_init, mu, ev_thresh, alpha, K);
[BG, FG, L_hat, S_hat, T_hat, t_hat] = ...
    ReProCS(M(:, t_train + 1 : end), ...
    L_init, mu, ev_thresh, alpha, K);
toc

L  = BG;
S = FG;



%VidName = ['Lobby_AutoReProCS_InitAltProj_alpha', num2str(alpha), ...
%    '_rank', num2str(rank_init)];
%DisplayVideo(M(:, t_train + 1 : end), FG, BG, T_hat, imSize, VidName);

