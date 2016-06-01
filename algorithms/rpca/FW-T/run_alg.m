% FW-T: SPCP solved by Frank-Wolfe method (Mu et al. 2014)
% process_video('RPCA', 'FW-T', 'dataset/demo.avi', 'output/demo_FW-T.avi');
[m,n] = size(M);
D = M/norm(M,'fro'); % imagesc(M); imagesc(D);

% parameter tuning
rho = 1; % rho = 0.5;  % sampling ratio
Omega = rand(m,n) <= rho; % support of observation imagesc(Omega);
obs = Omega.*D; % measurements imagesc(obs);

% this is parameter to control noise level
% the smaller the noise, the smaller is delta
delta = 0.01;

lambda_1 = delta*rho*norm(obs,'fro');
lambda_2 = delta*sqrt(rho)*norm(obs,'fro')/sqrt(max(m,n));

par.M = D;
par.lambda_1 = lambda_1;
par.lambda_2 = lambda_2;
par.iter = 1000;
par.display = 1;
par.rho = rho;
par.epsilon = 10^-3; % stopping criterion
par.method = 'exact'; % 'exact' or 'power'
par.Omega = Omega; % ones(m,n)
par.compare = 0; % make comparison or not

output = FW_T(par); % main function
%output = fista(par); % fista function
%output = ista(par); % ista function

L = output.L;
S = output.S;