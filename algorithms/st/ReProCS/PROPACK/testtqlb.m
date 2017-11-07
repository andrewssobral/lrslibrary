% Script for comparing speed an accuracy of original TQLB, optimized TQLB
% and builtin EIG command.

% Rasmus Munk Larsen, DAIMI, 1998.

n=1000;

% Use 2. order difference matrix as testproblem.
e = ones(n,1);
T = spdiags([-e 2*e -e], -1:1, n, n); 
true = 4*cos(pi/2*[n:-1:1]'./(n+1)).^2;
alpha = 2*ones(n,1);
beta = -ones(n,1);

fprintf('-----------------------------------------------------------------\n')
disp('Modified tqlb:')
fprintf('\n')
tic, flops(0)
[lambda,top,bot,err] = tqlb(alpha,beta);
fprintf('Elapsed time    = %f\n',toc);
fprintf('Number of flops    = %f\n',flops);
fprintf('Max rel. error  = %e\n',max(abs((lambda-true)./true)))


fprintf('-----------------------------------------------------------------\n')
disp('Original tqlb:')
fprintf('\n')
tic, flops(0);
[lambda2,top,bot,err2] = tqlb_orig(alpha,beta);
fprintf('Elapsed time    = %f\n',toc);
fprintf('Number of flops    = %f\n',flops);
fprintf('Max rel. error  = %e\n',max(abs((lambda2-true)./true)))


fprintf('-----------------------------------------------------------------\n')
disp('eig:')
fprintf('\n')
tic, flops(0);
lambda1 = eig(T);
lambda1 =sort(lambda1);
fprintf('Elapsed time    = %f\n',toc);
fprintf('Number of flops    = %f\n',flops);
fprintf('Max rel. error  = %e\n',max(abs((lambda1-true)./true)))

fprintf('-----------------------------------------------------------------\n')
disp('eig:')
fprintf('\n')
tic, flops(0);
lambda1 = eig(full(T));
lambda1 =sort(lambda1);
fprintf('Elapsed time    = %f\n',toc);
fprintf('Number of flops    = %f\n',flops);
fprintf('Max rel. error  = %e\n',max(abs((lambda1-true)./true)))


