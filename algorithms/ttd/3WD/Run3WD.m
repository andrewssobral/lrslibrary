% Three-Way Decomposition - Omar Oreifej 5/5/2012
clc ;
clear;
close all ;

datapath = 'C:\Research\Turbulence\DecompositionCode\Sequences\';

% input path
seqName = 'A' ;
imagePath = fullfile(datapath,[seqName, '\frames']) ;

% output path
destDir = fullfile(datapath,[seqName '\results3WD']);
if ~exist(destDir,'dir')
    mkdir(destDir) ;
end

% 3WD parameters
params.inner_tol = 1e-6 ;
params.maxIter = 1000 ;
params.mu = 1e-3 ;
params.tauc = 100*.04 ; % .5 for with Pi outside / 5 for with Pi inside. lambda = lambdac/sqrt(m) 1.1
params.lambdac = 1000;% 2000 good for Forb alone

[fileNames, numImages] = get_training_images(imagePath) ;

tmpim = imread(fileNames{1});
params.canonicalImageSize = size(tmpim);

% load object confidence file (can be obtained using any object detection
% method, including our proposed method)
load ([datapath '\' seqName '\ConfMap.mat']);
detThresh = .3; % .1 best for D, .3 best for A, .4 for C, .15 for B
[D, A, E,O, numIter] = prepare3WD(Probs,detThresh,fileNames, numImages, params,destDir);

saveResults(destDir, params.canonicalImageSize);

