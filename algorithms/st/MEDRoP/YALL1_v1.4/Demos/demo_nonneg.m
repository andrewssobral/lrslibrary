% solve the L1/L2 model with nonnegativity and compare with l1_ls
% edit the following line to access l1_ls
addpath(genpath('/Users/yin/Documents/matlab/l1_ls_matlab/'));
clear
n = 1000;
m = floor(.3*n);
k = floor(.2*m);
fprintf('\nSize [n,m,k] = [%i,%i,%i]\n',n,m,k);

% generate (A,b,xs) 
A = randn(m,n);
d = 1./sqrt(sum(A.^2));
A = A*sparse(1:n,1:n,d);
xs = zeros(n,1);
p = randperm(n); 
xs(p(1:k)) = 1; %rand(k,1);
b = A*xs;

% orthogonalize rows of A 
opts.nonorth = 0;
if ~opts.nonorth;
    [Q, R] = qr(A',0); 
    A = Q'; b = R'\b; 
end

% call YALL1 
opts.tol = 1e-6; 
opts.nonneg = 1;
disp('--- YALL1: BP');
tic; [x,Out] = yall1(A, b, opts); toc
fprintf('iter: %4i, rel_err = %6.2e\n',...
    Out.iter,norm(x-xs)/norm(xs));

disp('--- YALL1: L1/L1 (with wrong b)');
j = 3;
b(end-j:end) = -b(end-j:end);  % wrong b
opts.nu = .5;
tic; [x,Out] = yall1(A, b, opts); toc
fprintf('iter: %4i, rel_err = %6.2e\n',...
    Out.iter,norm(x-xs)/norm(xs));
b(end-j:end) = -b(end-j:end);  % restore b

sigma = .01;
rho = 1e-2;
tol = 1e-8;

disp('--- YALL1: L1L2');
b = b + sigma*randn(m,1);
opts.tol = tol; 
opts.nu = 0; 
opts.rho = rho; 
tic; [x1,Out] = yall1(A, b, opts); toc

% call l1_ls
disp('--- l1_ls');
if ~exist('l1_ls_nonneg','file'); error('Solver l1_ls_nonneg is not found.'); end
tic; x2 = l1_ls_nonneg(A, b, 2*rho, tol, 1); toc

fprintf('f_yall = %9.4f\n',norm(x1)+norm(A*x1-b)^2/2/rho);
fprintf('f_l1ls = %9.4f\n',norm(x2)+norm(A*x2-b)^2/2/rho);
fprintf('Rel_dist: %6.4e\n\n',norm(x1-x2)/norm(x1))

set(gca,'fontsize',18)
plot(1:n,x1,'bo',1:n,x2,'r.',1:n,xs,'k*');
legend('YALL1','l1-ls','Exact','location','best')
