
%[numr,numc] = size(M);
%I = randi([0 1],numr,numc);
I = ones(size(M));

mu = 1;
rho = 0.01;
toler = 1e-3;
maxiter = 100;
L = mc_logdet(M.*I,mu,rho,toler,maxiter);
S = M - L;
