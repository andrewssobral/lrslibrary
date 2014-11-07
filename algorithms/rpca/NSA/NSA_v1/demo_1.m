clc
clear all
addpath(genpath('./'));

tol = 5e-6; % optimality tolerance for stopping_type 1 
DB=20; % SNR of the video
frame_array=[35,100,125]; % frames that will be shown at the end 

seed = 602;
randn('state',seed); rand('state',seed);
if exist('Hall_airport_1000_1496_497_144_176_gray.mat')==0
    error('Video file cannot be found. Please download it from: http://www2.ie.psu.edu/aybat/Hall_airport_1000_1496_497_144_176_gray.mat')
else
    load Hall_airport_1000_1496_497_144_176_gray.mat;
end
D = images(:,1:201);
n1=144*176; n2=201;
stdev = norm(D,'fro')/(sqrt(144*176*201)*10^(DB/20));
D = D+stdev*randn(144*176,201);
[X,S]=nsa(D,stdev,tol,1);
figure
plot_data(frame_array,D,X,S,144,176)