%
% A variational approach to SPCP (Aravkin et al. 2014)
%
% RPCA | SPCP-sum-SPG | Stable PCP-sum solved by Spectral Projected Gradient (Aravkin et al. 2014)
% process_video('RPCA', 'flip-SPCP-sum-SPG', 'dataset/demo.avi', 'output/demo_flip-SPCP-sum-SPG.avi');

alg_path_aux = fullfile(lrs_conf.rpca_path,'SPGL1');
addpath(genpath(alg_path_aux));

nFrames     = size(M,2);
lambda      = 1/sqrt(max(size(M,1),size(M,2)));
L0          = repmat(median(M,2), 1, nFrames);
S0          = M - L0;
epsilon     = 5e-3*norm(M,'fro'); % tolerance for fidelity to data

% Flip-Flop version of SPCP-sum solved by Spectral Projected Gradient
opts = struct('sum',true,'L0',L0,'S0',S0,'max',false,...
  'tau0',3e5,'SPGL1_tol',1e-1,'tol',1e-3);
[L,S] = solver_RPCA_SPGL1(M,lambda,epsilon,[],opts);

rmpath(genpath(alg_path_aux));