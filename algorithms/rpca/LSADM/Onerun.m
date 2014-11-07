%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test ALM in the paper ''Fast Alternating Linearization Methods for
%       Minimizing the Sum of Two Convex Functions'', Donald Goldfarb,
%       Shiqian Ma and Katya Scheinberg, Tech. Report, Columbia University,
%       2009 - 2010. 
%
% Author: Shiqian Ma
% Date  : Apr. 20, 2010 
% IEOR, Columbia University, Copyright (2010)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

randn('state',0); rand('state',0);
% addpath('C:\Mywork\Optimization\work\code\ADM\Sparse_LowRank_Mat\Matdata');

% get data 
tic;
dataformat = 'surveillance-video-Hall';
opts = getdata(dataformat); 
time_getdata = toc;
fprintf('%f seconds to get data ! \n', time_getdata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Call ALM to solve the problem
tic; out_ALM = ALM_SADAL_smoothed(opts.D,opts); time_ALM = toc;
fprintf('*******************************************************************\n');
%%%%%%%%% print stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('*******************************************************************\n');
% fprintf('ALM  : iter: %d, relX: %3.2e, relY: %3.2e, StopCrit: %3.2e, time: %f\n', ...
%     out_ALM.iter, out_ALM.relX, out_ALM.relY, out_ALM.StopCrit, time_ALM);

% plot the figures

j = 10; 
subplot(3,3,1); imshow(reshape(opts.D(:,j),144,176),[]); 
subplot(3,3,2); imshow(reshape(out_ALM.X(:,j),144,176),[]); 
subplot(3,3,3); imshow(reshape(out_ALM.Y(:,j),144,176),[]); 

j = 80;
subplot(3,3,4); imshow(reshape(opts.D(:,j),144,176),[]); 
subplot(3,3,5); imshow(reshape(out_ALM.X(:,j),144,176),[]); 
subplot(3,3,6); imshow(reshape(out_ALM.Y(:,j),144,176),[]); 

j = 150;
subplot(3,3,7); imshow(reshape(opts.D(:,j),144,176),[]); 
subplot(3,3,8); imshow(reshape(out_ALM.X(:,j),144,176),[]); 
subplot(3,3,9); imshow(reshape(out_ALM.Y(:,j),144,176),[]); 