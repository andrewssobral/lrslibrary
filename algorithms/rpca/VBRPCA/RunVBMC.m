%% An example run of Variational Matrix Low-Rank Matrix Completion
clear
clc

%% Create a low-rank matrix
% Dimensions
m = 300;
n = 300;
r = 10; % rank

A = randn(m,r);
B = randn(r,n);
X_true = A*B; % low rank 

% random sparse sampling
df = r*(m+n-r); % Degrees of freedom
p = 0.2; % percentage measurements
pmn = p*m*n;
osdf = pmn / df; % Oversampling degrees of freedom

Omega = randperm(m*n); Omega = Omega(1:pmn);
P = zeros(m,n);  P(Omega) = 1;

% Observations
Y = P.*X_true;

% Add noise?
sigma = 0; % no noise
sigma = 1e-3; % some noise
% sigma = sqrt(1e-3); % more noise

Y = Y + sigma*randn(size(Y));

SNR = 10*log10( sum(sum((P.*X_true).^2)) / sum(sum((Y-P.*X_true).^2)));

%% Run VBMC
% all options are *optional*, everything will be set automatically
% you can modify these options to get better performance
options.verbose = 1;
options.MAXITER = 100; 
options.DIMRED = 1; % Reduce dimensionality during iterations?
% you can also set the threshold to reduce dimensions
% options.DIMRED_THR = 1e3;
%Estimate noise variance? (beta is inverse noise variance)
options.UPDATE_BETA = 1;
% options.UPDATE_BETA_START = 1;% iteration number to start estimating noise variance

% If the noise inv. variance is not to be estimated, set
% options.UPDATE_BETA = 0; % and set beta using
% options.beta = 1e3; 
% Manually tuning this parameter can give significantly better results. 

% options.initial_rank = 'auto'; % This sets to the maximum possible rank
options.initial_rank = 50; % or we can set a value. 

% Provide original matrix for simulations
options.X_true = X_true; 



fprintf('Running VB matrix completion...\n');
tic
[X_hat, A_hat, B_hat] = VBMC(P, Y, options);
t_total = toc;

% Show results
err = norm( X_hat - X_true, 'fro' ) / norm( X_true, 'fro');
r_hat = rank(X_hat);
fprintf('\n\n-------------------------------------------------------------------\n')
fprintf('Dimensions (%d, %d), true rank = %d, measurements = %g %%, SNR = %g\n',m, n, r, p*100, SNR);
fprintf('Rel. error = %g, estimated rank = %d, time = %g\n', err, r_hat, t_total);
fprintf('-------------------------------------------------------------------\n')


