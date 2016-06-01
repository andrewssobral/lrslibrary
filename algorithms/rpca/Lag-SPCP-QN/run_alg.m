%
% A variational approach to SPCP (Aravkin et al. 2014)
%
% RPCA | Lag-SPCP-QN  | Lagrangian SPCP solved by Quasi-Newton (Aravkin et al. 2014)
% process_video('RPCA', 'Lag-SPCP-QN', 'dataset/demo.avi', 'output/demo_Lag-SPCP-QN.avi');

alg_path_aux = fullfile(lrs_conf.rpca_path,'SPGL1');
addpath(genpath(alg_path_aux));

nFrames     = size(M,2);
lambda      = 1/sqrt(max(size(M,1),size(M,2)));
L0          = repmat(median(M,2), 1, nFrames);
S0          = M - L0;
epsilon     = 5e-3*norm(M,'fro'); % tolerance for fidelity to data

% Lagrangian SPCP solved by Quasi-Newton
opts    = struct('sum',false,'L0',L0,'S0',S0,'max',false,'tol',1e-3,'quasiNewton',true);
lambdaL = 0.25; lambdaS = 0.01;
[L,S] = solver_RPCA_Lagrangian(M,lambdaL,lambdaS,[],opts);

rmpath(genpath(alg_path_aux));