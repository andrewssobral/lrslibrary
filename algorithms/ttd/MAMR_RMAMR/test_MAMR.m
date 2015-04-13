%test_MAMR.m tests the condition of binary mask in Omega.
%input: car_data.mat(frames) car_mask.mat(obtained from optical flow(binary))
% Copyright: Xinchen YE, Tianjin University, 2014
clear;clc;addpath('PROPACK');
load car_data.mat 
load car_mask.mat
% show_2dvideo(D,240,320)
% show_2dvideo(confidence_mask,240,320)
noise_v = 10; %noise variance added on frames in D
lambda = 1; 
[m, n]=size(D);
D1 = D + noise_v * randn(size(D));  % show_2dvideo(D1,240,320)
D2 = D1 .* confidence_mask; % show_2dvideo(D2,240,320)
Omega = find(D2 ~= 0);
[I, J] = ind2sub([m n],Omega);  
%%
tic
   %[L,S,iter1]=core_MAMR(D, lambda, I, J);
   [L,S,Z,iter1] = core_RMAMR(D,10, I, J);
toc
% psnr = calcpsnr(reshape(F(:,12),288,352), double(monitor_back_gray))
%%
% show_2dvideo(L,240,320) show_2dvideo(S,240,320) show_2dvideo(Z,240,320)
for i = 1:24
    figure(101);imshow([reshape(D(:,i),240,320) reshape(L(:,i),240,320)],[]); pause(0.1);
end