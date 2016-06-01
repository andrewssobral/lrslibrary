%Test script to run SVP for matrix completion using uniformly sampled random matrices
%Written by: Prateek Jain (pjain@cs.utexas.edu) and Raghu Meka (meka@cs.utexas.edu)
%Last updated: October, 2009

m = 5000;
n = 2000;
k = 10; %rank of the optimal matrix
p = 0.2; %fraction of known entries
params.tol=1e-3;
params.vtol=1e-4;
params.mxitr=100;
params.verbosity=1;

M = randn(m,k)*randn(k,n); %Generate the underlying matrix
random_permutation=randperm(m*n);
Omega=random_permutation(1:p*m*n)'; %Generate samples with uniform density p
Omega=sort(Omega);

fprintf('Running SVP with m=%d, n=%d, p=%f...\n',m,n,p);
t=cputime;
[U,S,V,num_iter] = svp(Omega, M(Omega), m, n, k, params); %Run SVP
t = cputime-t;
fprintf('SVP finished with RMSE=%f in time t=%f\n',norm(M-U*diag(S)*V','fro')/sqrt(length(Omega)),t);
