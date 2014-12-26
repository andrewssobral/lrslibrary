%demo for the dinosaur sequence

clear M;
clear saturation;
clear occlusion;

load M_yaleb10.mat
load saturation_yaleb10.mat
load occlusion_yaleb10.mat

disp('-------------------------------Low-Rank Matrix Factorization Under L1-Norm------------------------------------')
disp('------------------------Rank-3 Factorization for Extended Yale Face B10 Sequence------------------------------')
fprintf('\n');

row = 192;
col = 168;

%indicator matrix
W = saturation.*occlusion;

tic;
[M_est U_est V_est L1_error] = RobustApproximation_M_UV_TraceNormReg(M,W,3,1e-3,1.2,1,1);

time = toc;
%draw the results for the 46-th frame
index = 46;
figure;
subplot(1,3,1);imshow(reshape(M(index,:),row,col)); title(['Input-',num2str(index),'th frame']);
vim(:,:,2) = reshape(occlusion(index,:),row,col);
vim(:,:,1) = ones(row,col);
vim(:,:,3) = reshape(saturation(index,:),row,col);
subplot(1,3,2); imshow(vim); title('Mask');
subplot(1,3,3); imshow(reshape(M_est(index,:),row,col)); title('Recovery');

fprintf('\n');
disp('-----------------------Performance Statistics--------------------- ');
%show average L1_error
disp(['Average L1-Norm Error (8-bit depth): ', num2str(255*L1_error/sum(sum(W)))]);

%show running time
disp(['Total Runing Time: ', num2str(time), 'seconds']);
