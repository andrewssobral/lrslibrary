function demo2
    clc
    clear
    close all

    addpath('./util/build-group/');
    addpath('./util');

    regularizer = 0.1;
    param.superpixel.slicParam = [10, regularizer; 20, regularizer;  40,regularizer ; 80,regularizer; ];

    %  param.superpixel.slicParam = [5,regularizer;  10, regularizer; 20, regularizer;  40,regularizer ; ];

    param.maxIter =100;
    param.tol = 1e-4;

    param.input.dataPath = '/Users/Jia/research/dataset/wallflower/';
    param.output.resultPath = '/tmp/result-wallflower';
    param.input.dirString = 'b*.bmp';

    param.input.dataset = 'MovedObject'; 

    param.saveResult = 0;
    param.showResult = 1;

    param.output.threshold  = 0.1;
    param.lambda =10;   

    param.rank = 5; % subspace dimension
    param.eta = 0.01;  % stepsize for subspace updating
%     param.sampleSize = 200; % number of sample frames to learn the background
%     param.trainSize = 5000; % train with all available images
    param.startIndex = 1;
    param.randomStart = true;
  
        
    param     
     
    gosus(param);  

 

end



 