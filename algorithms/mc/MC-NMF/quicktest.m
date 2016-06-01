clear

% problem specification
m = 400; n = 500; r = 10;
t = 1; esr = ceil(t*r); 
fprintf('[m n r esr] = [%i, %i, %i, %i]\n',m,n,r,esr);

% data
% Xo = rand(m,r); Yo = rand(r,n);
Xo = abs(randn(m,r)); Yo = abs(randn(r,n));

d = ones(r,1).*(1:r)';
M = Xo*spdiags(d,0,r,r)*Yo;

% set solver options
opts.tol = 1e-5;
opts.maxit = 500;
opts.print = 1;

% call solver

sr = 0.5;
Omega = randsample(m*n, sr*m*n);
A = M(Omega);

%exact rank-estimate case
t0 = tic;
[X,Y,Out] = mc_nmf(A,Omega,esr,m,n,opts);
time = toc(t0);
relerr = norm(X*Y-M,'fro')/norm(M,'fro');
fprintf('RelErr = %5.2e, %4.2e Sec.\n',relerr,time);
%%
[numr,numc] = size(M);
I = randi([0 1],numr,numc);
Omega = find(I);
A = M(Omega);

%%
subplot(1,2,1),imagesc(M);
subplot(1,2,2),imagesc(X*Y);
%%

%over-estimate case
esr = ceil(1.5*r);
t0 = tic;
fprintf('[m n r esr] = [%i, %i, %i, %i]\n',m,n,r,esr);
[X,Y,Out] = mc_nmf(A,Omega,esr,m,n,opts);
time = toc(t0);
relerr = norm(X*Y-M,'fro')/norm(M,'fro');
fprintf('RelErr = %5.2e, %4.2e Sec.\n',relerr,time);