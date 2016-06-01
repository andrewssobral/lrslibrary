% This is the implementation of nonconvex RPCA method in ICDM 2015 paper
% 'Robust PCA via Nonconvex Rank Approximation'
% Zhao Kang, August 2015. Questions? Zhao.Kang@siu.edu;

[m,n] = size(M);
muzero = .9; % the only tuning parameter
lambda = 1/sqrt(max(size(M))); % default lambda
%lambda = 1e-3;  % model parameter
type = 21;   %different modeling of Sparsity.
rate = 1.1;   %update rate of \mu
gamma = 1e-2;     %gamma parameter in the rank approximation
tol = 1e-6;  % stopping criterion

%initializations
S = zeros(m,n);
Y = zeros(m,n);
L = M;
sig = zeros(min(m,n),1); % for DC
mu = muzero;

for ii=1:500
  D = M-S-Y/mu;
  [L,sig] = DC(D,mu/2,sig,gamma);
  [S] = errorsol(Y,M,L,lambda,mu,type);
  Y = Y + mu*(L-M+S);
  mu = mu * rate;
  
  sigma = norm(M-S-L,'fro');
  RRE = sigma/norm(M,'fro');
  
  if RRE < tol
    break;
  end  
end
rk = rank(L);
disp(['  relative err ',num2str(RRE),' rank ',num2str(rk)]);
