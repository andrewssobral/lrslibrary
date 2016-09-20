
%[numr,numc] = size(M);
%Idx = randi([0 1],numr,numc);
Idx = ones(size(M));

mu = 1;
rho = 0.01;
toler = 1e-3;
maxiter = 100;
L = mc_logdet(M.*Idx,mu,rho,toler,maxiter);
S = M - L;
