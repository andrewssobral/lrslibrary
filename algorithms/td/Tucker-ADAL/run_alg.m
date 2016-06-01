% TD | Tucker-ADAL | Tucker Decomposition solved by ADAL (Goldfarb and Qin, 2013)
% process_video('TD', 'Tucker-ADAL', 'dataset/demo.avi', 'output/demo_Tucker-ADAL.avi');

alg_path_aux = fullfile(lrs_conf.td_path,'RLRT');
addpath(genpath(alg_path_aux));

pdata.T = T;
pdata.X = T;
N = ndims(pdata.T);
r = 1/sqrt(max(size(pdata.T)));
params.E0 = tenzeros(size(pdata.T));
params.X0 = tenzeros(size(pdata.T));
params.V0 = cell(1, N);
for i = 1:N
  params.V0{i} = tenzeros(size(pdata.T));
end
params.mu0 = 1/(N+1);
params.mode = N;
params.IsTC = false; % is tensor completion
params.rRatio = 1/4;
params.opt_tol = 1e-3;
params.eta = 1/(N+1);
params.max_iter = 1000;
params.mu1fac = 10;
params.mu1 = params.mu1fac*std(T(:));
params.mu2 = params.mu1;
params.mu_min = 1e-4;
params.mu_max = 1e2;
params.lambdaS = 1;
params.lambda = params.lambdaS*r*params.rRatio;
params.verbose = 1;
params.use_cont = true;
params.k = [size(T,1) size(T,2) 1];
%%%%%%%%%% for PROPACK %%%%%%%%%%%%
% declare global var 'sv'
global sv;
global tmode;
global use_propack;
global curr_mu;
sv =  ceil(min(size(pdata.T)) * 0.1) * ones( 1, N );
use_propack = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
results = tensor_tucker_adal_ncx(pdata, params);
L = double(results.X);
S = double(results.E);
clear sv tmode use_propack curr_mu;

rmpath(genpath(alg_path_aux));
