% solve the L1/L2 model and compare with l1_ls
% edit the following line to access l1_ls
addpath(genpath('/Users/yin/Documents/matlab/l1_ls_matlab/'));
clear
n = 400;
m = floor(.4*n);
k = floor(.2*m);
fprintf('\nSize [n,m,k] = [%i,%i,%i]\n',n,m,k);

% generate (A,b,xs) 
sigma = 0.5;
A = randn(m,n);
xs = zeros(n,1); 
p = randperm(n); 
xs(p(1:k)) = 2*(rand(k,1) > 0.5) - 1;
b = A*xs + sigma*randn(m,1); 

% orthogonalize rows of A
opts.nonorth = randn > 0;
fprintf('nonorth = %i\n',opts.nonorth);
if opts.nonorth;
    d = 1./sqrt(sum(A.^2,2));
    A = sparse(1:m,1:m,d)*A;
    b = b.*d;
else
    [Q, R] = qr(A',0);
    A = Q'; b = R' \ b;
end
rho = 1e-3;

% call YALL1 
disp('--- YALL1 ---');
opts.tol = 1e-3; opts.rho = rho;
tic; [x1,Out] = yall1(A, b, opts); toc

% call l1_ls
disp('--- l1_ls ---');
if ~exist('l1_ls','file'); error('Solver l1_ls is not found.'); end
tic; x2 = l1_ls(A, b, 2*rho, opts.tol, 1); toc

dx12 = norm(x1-x2)/norm(x1);
fprintf('Diff. of the two: %6.4e\n\n',dx12)

set(gca,'fontsize',18)
plot(1:n,x1,'bo',1:n,x2,'r.',1:n,xs,'k*');
legend('YALL1','l1-ls','Exact')
