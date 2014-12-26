clc
clear
close all

%   I = imread('~/research/dataset/stat-background-substraction/WaterSurface/WaterSurface1632.bmp');
 
%    I = imread('~/research/dataset/stat-background-substraction/Campus/trees1792.bmp');

 
I = imread('/home/jxu/research/dataset/wallflower-background-substraction/WavingTrees/b00247.bmp');


 
% I = imread('./41006.jpg');
 % I = imread('./b00247.bmp');
 
 

%  segment = vl_slic(I,10,2);
 imlab = vl_xyz2lab(vl_rgb2xyz(single(I))) ;
 imlab = single(imlab);
  tic
segments = vl_slic(imlab, 5, 0.1) ;
 toc
 
 



  tic
segments = vl_slic(imlab, 20, 0.1) ;
 toc
 
 
tic
segments = vl_slic(imlab, 20, 0.01) ;

 toc