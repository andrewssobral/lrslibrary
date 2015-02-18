% A demo showing how to use DRMF for video background extraction.
% Results are compared to RPCA by Candes et al.
%
% You should be able to run this demo directly in Matlab.
% 
% author: Liang Xiong (lxiong@cs.cmu.edu)

addpath ./propack
addpath ./yima_rpca
addpath ./lib

load fgbg_restaurant200;
X = imgs; clear imgs;
[m, n] = size(X);

%% RPCA

disp('RPCA');

lambda = 1/sqrt(max(size(X)));
ts = tic;
[L_rpca S_rpca] = inexact_alm_rpca(X, lambda);
t_rpca = toc(ts);
% ** Our later study shows that debiasing should not be applied 
%    in the outlier detection situations, unlike the original paper.

% guess the rank based on variance preservation
sv = svdex(L_rpca);
rk = EffRank(sv, 0.999);

%% SVD

disp('SVD');

ts = tic;
[U_svd S_svd V_svd] = svdex(X, rk);
t_svd = toc(ts);

L_svd = bsxfun(@times, U_svd, S_svd')*V_svd';
S_svd = X - L_svd;

%% DRMF

disp('DRMF');

% initialize DRMF to get to a better local minimum.
options.init = inexact_alm_rpca(X, lambda, 1e-5, 10);

ts = tic;
[L_drmf S_drmf] = DRMF(X, rk, 0.1, options);
t_drmf = toc(ts);

% cache the results
save cache_drmf_demo

%% report the running time

fprintf('Running time: \nRPCA = %0.1f sec\nSVD  = %0.1f sec\nDRMF = %0.1f sec\n',...
    t_rpca, t_svd, t_drmf);

%% compare different methods on the sequence

figure;
for ind = 1:size(X,1)
    ShowImages([X(ind,:)' L_svd(ind,:)' L_rpca(ind,:)' L_drmf(ind,:)'],imsize,4);
    title(strrep('Original frame | SVD background | RPCA background | DRMF background',...
        '|', '            |            '));
    drawnow;    pause(0.1);
end

%% compare different methods on selected frames

% randomly select some frames to compare
idx = ceil(rand(5,1)*size(X,1));

orig = X(idx,:);
bg_svd = L_svd(idx,:);
bg_rpca = L_rpca(idx,:);
bg_drmf = L_drmf(idx,:);

figure;
for ind = 1:100
    ShowImages([orig' bg_svd'], imsize, length(idx));
    ylabel(strrep('Background | Original Frame','|','        |        '));
    title('Method: SVD');    pause(0.5);
    
    ShowImages([orig' bg_rpca'], imsize, length(idx));
    ylabel(strrep('Background | Original Frame','|','        |        '));
    title('Method: RPCA');   pause(0.5);
    
    ShowImages([orig' bg_drmf'], imsize, length(idx));
    ylabel(strrep('Background | Original Frame','|','        |        '));
    title('Method: DRMF');   pause(0.5);
end
