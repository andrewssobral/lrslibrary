% restoredefaultpath;
function lrs_setup()
disp('Running setup...');

[lrs_dir,~,~] = fileparts(mfilename('fullpath'));
lrs_conf.lrs_dir = lrs_dir;
disp(['LRSLibrary path: ' lrs_dir]);

lrs_conf.lrr_path  = fullfile(lrs_dir,'algorithms','lrr');
lrs_conf.nmf_path  = fullfile(lrs_dir,'algorithms','nmf');
lrs_conf.ntf_path  = fullfile(lrs_dir,'algorithms','ntf');
lrs_conf.rpca_path = fullfile(lrs_dir,'algorithms','rpca');
lrs_conf.td_path   = fullfile(lrs_dir,'algorithms','td');

lib_path = fullfile(lrs_dir,'libs');
lrs_conf.pro_path = fullfile(lib_path,'PROPACK2');
lrs_conf.svd_path = fullfile(lib_path,'SVD');
lrs_conf.vl_path = fullfile(lib_path,'vlfeat-0.9.19','toolbox','vl_setup');
lrs_conf.tfocs_path = fullfile(lib_path,'TFOCS-1.3.1');
lrs_conf.cvx_path = fullfile(lib_path,'cvx');
lrs_conf.tensor_toolbox_path = fullfile(lib_path,'tensor_toolbox_2.5');
lrs_conf.mtt_path = fullfile(lib_path,'mtt');
lrs_conf.nway_path = fullfile(lib_path,'nway331');
lrs_conf.poblano_toolbox_path = fullfile(lib_path,'poblano_toolbox_1.1');
lrs_conf.lightspeed_path = fullfile(lib_path,'lightspeed');
lrs_conf.tensorlab_path = fullfile(lib_path,'tensorlab');
lrs_conf.lib_path = lib_path;

save(fullfile(lrs_dir,'lrs_conf.mat'),'lrs_conf');

disp('Updating PATH!');
addpath(lrs_conf.lrs_dir); % add LRSLibrary
addpath(lrs_conf.pro_path); % add PROPACK
addpath(lrs_conf.svd_path); % add SVD
addpath(lrs_conf.tfocs_path); % add TFOCS
addpath(lrs_conf.tensor_toolbox_path); % add Tensor-Toolbox
addpath(lrs_conf.mtt_path); % add MTT
addpath(lrs_conf.nway_path); % add NWay
addpath(lrs_conf.poblano_toolbox_path); % add Poblano-Toolbox
addpath(lrs_conf.lightspeed_path); % add Lightspeed
addpath(lrs_conf.tensorlab_path); % add TensorLab

disp('Running VLFeat setup!');
run(lrs_conf.vl_path); % add VLFeat

if(~isemptydir(lrs_conf.cvx_path))
  addpath(lrs_conf.cvx_path); % add CVX
  disp('Running CVX setup!');
  cvx_setup;
end

disp('Setup finished!');
end