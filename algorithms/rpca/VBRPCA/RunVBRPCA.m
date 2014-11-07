%% An example run of VBRPCA
clear
clc

%% Create a low-rank + sparse matrix
% Dimensions
m = 500;
n = 500;
r = 25; % rank of the low-rank component

A = randn(m,r);
B = randn(r,n);
X_true = A*B; % low rank component

SparseRatio = 0.05;
E_true = zeros(size(X_true));
Eomega = randsample(m*n, round(m*n*SparseRatio));
E_true(Eomega) = -10+20*rand(length(Eomega),1); % sparse component

% Observation
Y = X_true + E_true;

% Add noise?
sigma = 0; % no noise
% sigma = 1e-3; % some noise
% sigma = sqrt(1e-3); % more noise

Y = Y + sigma*randn(size(Y));

SNR = 10*log10( sum(sum((X_true+E_true).^2)) / sum(sum((Y-X_true-E_true).^2)));


%% Run VRBPCA
% all options are *optional*, everything will be set automatically
% you can modify these options to get better performance
options.verbose = 1;
options.initial_rank = 'auto'; % This sets to the maximum possible rank
% options.initial_rank = 300; % or we can use a value. 
options.X_true = X_true;
options.E_true = E_true;
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


fprintf('Running VBRPCA...\n')
tic
[X_hat, A_hat, B_hat, E_hat] = VBRPCA(Y,options);
t_total = toc;

% Results
Xerr = norm( X_hat - X_true, 'fro' ) / norm( X_true, 'fro');
Eerr = norm( E_hat - E_true, 'fro' ) / norm( E_true, 'fro');

r_hat = rank(X_hat);
% the model does not allow for "exact" sparsity, but only "soft" sparsity.
% That is, most coefficients are almost zero (very small values), but not 
% exactly zero. If the sparse coefficients are needed, one can use
% thresholding such as
nonzero_E = length(find(abs(E_hat)>1e-7));

% Show results
fprintf('\n\n-------------------------------------------------------------------\n')
fprintf('Dimensions (%d, %d), true rank = %d, nonzero E elements = %d, SNR = %g\n',m, n, r, length(Eomega), SNR);
fprintf('low rank comp. error = %g, sparse comp. error %g, est. rank = %d,\nnonzero E elements = %d, running time = %g\n', Xerr, Eerr, r_hat, nonzero_E, t_total);
fprintf('-------------------------------------------------------------------\n')

