%clear;

if ~exist('maskmult')
  fprintf('Mex functions are not available\n');
  compile_mex;
end

m = 5000;
n = 5000;
k = 5;
eps = 0.01; % fraction of observed entries
noise = 0.0;
fprintf('Generating matrix m=%d, n=%d, k=%d, observed percentage = %.2f%%, noise = %f\n', m, n, k, eps*100, noise);

%%---------- Generate data

% generate random low rank matrix A
test_case = 1;
switch test_case
case 1
  % Test case 1: random low rank matrix with chosen singular values
  [U0 junk] = qr(randn(m,k),0); [V0 junk] = qr(randn(n,k),0);
  s1 = 1000; step=1000; S0 = diag([s1:step:s1+(k-1)*step]*1);
  A0 = U0*S0*V0'; 
case 2
  % Test case 2: random low rank matrix
  U0 = rand(m,k); V0 = rand(n,k); A0 = U0*V0';
end

% observed entries at random positions
K = rand(m,n)<=eps;
N = (K==0);
% add noise
A = A0 + noise*(K.*randn(m,n));
% get Train and Test data
Train = sparse(K.*A);
Test = N.*A;
[I J VALS] = find(Train);
I = int32(I);
J = int32(J);

%%---------- Run the methods ----------

% array to save history data
meth_name = {};
meth_obj = {};
meth_time = {};

% run ScGrassMC
tol = 1.e-6;
tol_reschg = tol/10;
maxit = 100;
fprintf('ScGrassMC\n');
bt = tic;
% type help ScGrassMC for options for ScGrassMC
[U S V hist] = ScGrassMC(K, VALS, k,...
                  'tol', tol,...
                  'maxit',maxit,...
                  'grad_type','scaled',...
                  'beta_type','P-R',...
                  'sigma_type','approx',... %'tol_reschg', tol_reschg,...
                  'verbose', 1);
meth_name{end+1} = 'ScGrassMC';
meth_obj{end+1} = hist.rmse;
meth_time{end+1} = hist.time;
toc(bt);
fprintf('Test RMSE: %f\n', norm(N.*(U*S*V')-Test,'fro')/sqrt(nnz(N)));

% run other methods here
% ...

% plot results
fname = ''; % filename to save image, empty string to ignore saving
singvals = ''; % print singular values to put into graph's title
for i=1:k-1
  singvals = sprintf('%s%.0f,', singvals, S0(i,i));
end
singvals = sprintf('%s%.0f', singvals, S0(k,k));
title = sprintf('%dx%d - Rank %d - %.1f%% observed entries\nSingular values [%s]',m,n,k,eps*100,singvals);
plot_error_time(fname, meth_obj, meth_time, meth_name, title, 'log');
