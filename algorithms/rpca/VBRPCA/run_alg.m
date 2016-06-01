% Variational Bayesian RPCA (Babacan et al. 2011)
% process_video('RPCA', 'VBRPCA', 'dataset/demo.avi', 'output/demo_VBRPCA.avi');

% all options are *optional*, everything will be set automatically
% you can modify these options to get better performance
options.verbose = 1;
options.initial_rank = 'auto'; % This sets to the maximum possible rank
% options.initial_rank = 300; % or we can use a value.
%options.X_true = X_true;
%options.E_true = E_true;
options.inf_flag = 2; % inference flag for the sparse component
% 1 for standard VB, 2 for MacKay. MacKay generally converges faster.
options.MAXITER = 200;
%Estimate noise variance? (beta is inverse noise variance)
options.UPDATE_BETA = 1;
% If the noise inv. variance is not to be estimated, set
% options.UPDATE_BETA = 0; % and set beta using
% options.beta = 1e3;
% Select the optimization mode:
% 'VB': fully Bayesian inference (default)
% 'VB_app': fully Bayesian with covariance approximation
% 'MAP': maximum a posteriori (covariance is set to 0)
options.mode = 'VB';
% For large scale problems, set to 'VB_app'.
% options.mode = 'VB_app';
[L,~,~,S] = VBRPCA(M, options);
