% restoredefaultpath;
function lrs_setup()
disp('Running setup...');

[lrs_dir,~,~] = fileparts(mfilename('fullpath'));
lrs_conf.lrs_dir = lrs_dir;
disp(['LRSLibrary path: ' lrs_dir]);

lrs_conf.utils_path = fullfile(lrs_dir,'utils');
lrs_conf.gui_path  = fullfile(lrs_dir,'gui');

lrs_conf.lrr_path  = fullfile(lrs_dir,'algorithms','lrr');
lrs_conf.mc_path   = fullfile(lrs_dir,'algorithms','mc');
lrs_conf.nmf_path  = fullfile(lrs_dir,'algorithms','nmf');
lrs_conf.ntf_path  = fullfile(lrs_dir,'algorithms','ntf');
lrs_conf.rpca_path = fullfile(lrs_dir,'algorithms','rpca');
lrs_conf.st_path   = fullfile(lrs_dir,'algorithms','st');
lrs_conf.td_path   = fullfile(lrs_dir,'algorithms','td');
lrs_conf.ttd_path  = fullfile(lrs_dir,'algorithms','ttd');

lib_path = fullfile(lrs_dir,'libs');
lrs_conf.pro_path = fullfile(lib_path,'PROPACK');
lrs_conf.svd_path = fullfile(lib_path,'SVD');
lrs_conf.vl_path = fullfile(lib_path,'vlfeat','toolbox','vl_setup');
lrs_conf.tfocs_path = fullfile(lib_path,'TFOCS');
lrs_conf.yall1_path = fullfile(lib_path,'YALL1');
lrs_conf.cvx_path = fullfile(lib_path,'cvx');
lrs_conf.tensor_toolbox_path = fullfile(lib_path,'tensor_toolbox');
lrs_conf.mtt_path = fullfile(lib_path,'mtt');
lrs_conf.nway_path = fullfile(lib_path,'nway');
lrs_conf.poblano_toolbox_path = fullfile(lib_path,'poblano_toolbox');
%lrs_conf.lightspeed_path = fullfile(lib_path,'lightspeed');
%lrs_conf.tensorlab_path = fullfile(lib_path,'tensorlab');
%lrs_conf.mmread_path = fullfile(lib_path,'mmread');
lrs_conf.manopt_path = fullfile(lib_path,'manopt');
lrs_conf.lib_path = lib_path;

save(fullfile(lrs_dir,'lrs_conf.mat'),'lrs_conf');

disp('Updating PATH!');
addpath(lrs_conf.utils_path); % add LRSLibrary utils
addpath(lrs_conf.gui_path); % add LRSLibrary gui

addpath(lrs_conf.lrs_dir); % add LRSLibrary
addpath(lrs_conf.lib_path); % add Library Path
addpath(lrs_conf.pro_path); % add PROPACK
addpath(lrs_conf.svd_path); % add SVD
addpath(lrs_conf.tfocs_path); % add TFOCS
addpath(lrs_conf.yall1_path); % add YALL1
addpath(lrs_conf.tensor_toolbox_path); % add Tensor-Toolbox
addpath(lrs_conf.mtt_path); % add MTT
addpath(lrs_conf.nway_path); % add NWay
addpath(lrs_conf.poblano_toolbox_path); % add Poblano-Toolbox
%addpath(lrs_conf.lightspeed_path); % add Lightspeed
%addpath(lrs_conf.tensorlab_path); % add TensorLab
%addpath(lrs_conf.mmread_path); % add MMRead
addpath(lrs_conf.manopt_path); % add Manopt and subfolders
addpath(genpath(fullfile(lrs_conf.manopt_path,'manopt')));

disp('Running VLFeat setup!');
run(lrs_conf.vl_path); % add VLFeat

if(~isemptydir(lrs_conf.cvx_path,3) && exist(fullfile(lrs_conf.cvx_path,'cvx_setup.m'),'file'))
  addpath(lrs_conf.cvx_path); % add CVX
  disp('Running CVX setup!');
  cvx_setup;
end

disp('Setup finished!');
end