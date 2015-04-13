function demo1
    clc
    clear
    close all

    addpath('./util/build-group/');
    addpath('./util');


    regularizer = 0.1;
    param.superpixel.slicParam = [10, regularizer; 20, regularizer;  40,regularizer ; 80,regularizer; ];

    %  param.superpixel.slicParam = [5,regularizer;  10, regularizer; 20, regularizer;  40,regularizer ; ];

    param.maxIter =200;
    param.tol = 1e-4;

    % data path
    param.input.dataPath = '/Users/Jia/research/dataset/wallflower/';
    param.output.resultPath = '/tmp/result-wallflower';
    param.input.dirString = 'b*.bmp';

    param.input.dataset = 'WavingTrees'; 

    param.saveResult = 0;
    param.showResult = 1;

    param.output.threshold  = 8e-6;
    param.lambda =0.33;    %   bestAccu = 0.9762, lambda = 0.33; threshold = 8e-6;
    
    param.rank = 5; % subspace dimension
    param.eta = 1e-3 ;   % stepsize for subspace updating
    param.sampleSize = 200; % number of sample frames to learn the background
    param.trainSize = 240; % sample from the fist 240 frames to learn the background
    param.startIndex = 1;   %  which frame to start
    param.randomStart = false;  

 
    param

    gosus(param);  


end



