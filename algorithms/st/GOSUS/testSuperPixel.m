% Test Slic
clc
clear
close all

%    I = imread('~/research/dataset/stat-background-substraction/WaterSurface/WaterSurface1632.bmp');
 %  I = imread('~/research/dataset/stat-background-substraction/Fountain/Fountain1184.bmp');
 
%    I = imread('~/research/dataset/stat-background-substraction/Campus/trees1792.bmp');
 
  I = imread('/home/jxu/research/dataset/wallflower-background-substraction/WavingTrees/b00247.bmp');
 
 
 
 
 figure;
 subplot(1, 4, 1);
 imshow(I);

 
 tic
%  segment = vl_slic(I,10,2);
 imlab = vl_xyz2lab(vl_rgb2xyz(single(I))) ;
 imlab = single(imlab);
segments = vl_slic(imlab, 60, 0.1) ;
 toc
 
%  segments = single(segments);
% I_sp = segImage(I,segments);
I_sp = segImageRegion(I,segments);

subplot(1,4,2)
max(I_sp(:))
min(I_sp(:))
imshow(I_sp);

imwrite(I_sp, 'image/sp_region1.jpg', 'jpg');

 
segments = vl_slic(imlab, 20, 0.1) ;

 
%  segments = single(segments);
% I_sp = segImage(I,segments);
I_sp = segImageRegion(I,segments);

subplot(1,4,3)
% max(I_sp(:))
% min(I_sp(:))
imshow(I_sp);

imwrite(I_sp, 'image/sp_region2.jpg', 'jpg');
 

segments = vl_slic(imlab, 50, 0.1) ;

 
 
%  segments = single(segments);
% I_sp = segImage(I,segments);
I_sp = segImageRegion(I,segments);

subplot(1,4,4)
max(I_sp(:))
min(I_sp(:))
imshow(I_sp);

imwrite(I_sp, 'image/sp_region3.jpg', 'jpg');

