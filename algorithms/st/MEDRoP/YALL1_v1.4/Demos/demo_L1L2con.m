% solve the L1/L2con model and compare with SPGL1
% edit the following line to access SPGL1

addpath(genpath('/Users/yin/Documents/matlab/spgl1-1.7/'));
clear;

n = 1000;
m = round(.2*n);
k = round(.25*m);
fprintf('\nSize [n,m,k] = [%i,%i,%i]\n',n,m,k);

% generate (A,b,xs) 
A = randn(m,n);
xs = zeros(n,1); 
p = randperm(n); 
xs(p(1:k)) = 4*(rand(k,1)>.5) - 2;  
% nonnegative signal
xs = abs(xs);

sigma = 0.001;
noise = sigma*randn(m,1);
b = A*xs + noise;
delta = norm(noise);
fprintf('sigma = %6.3e, delta = %6.3e\n\n',sigma,delta);

if all(xs >= 0); J = [1 0]; else J = 0; end
% orthogonalize rows of A
[Q, R] = qr(A',0);
A = Q'; b = R' \ b;

% call YALL1 
disp('--- YALL1 ---');
opts.tol = 1e-3; 
opts.delta = delta;
opts.print = 1; 
for j = J
    opts.nonneg = j; fprintf('nonneg set to %i\n',j);
    tic; [x1,Out] = yall1(A, b, opts); toc
    fprintf('||Ax-b|| = %6.3e: ||x||_1 = %6.3e\n',...
        norm(A*x1-b), norm(x1,1));
    rerr = norm(x1-xs)/norm(xs);
    fprintf('[nA,nAt]=[%i,%i]: Rel_err = %6.2e\n\n',...
        Out.cntA,Out.cntAt,rerr);
end

% call SPGL1
disp('--- SPGL1 ---');
if ~exist('spgSetParms','file'); error('Solver SPGL1 is not found.'); end
spg_opts = spgSetParms('verbosity',0,'optTol',opts.tol);
tic; [x2,r,g,info] = spgl1(A,b,0,delta,[],spg_opts); toc
fprintf('||Ax-b|| = %6.3e: ||x||_1 = %6.3e\n',...
    norm(A*x2-b), norm(x2,1));
rerr = norm(x2-xs)/norm(xs);
fprintf('[nA,nAt]=[%i,%i]:  Rel_err = %6.2e\n\n',...
    info.nProdA,info.nProdAt,rerr);

dx12 = norm(x1-x2)/norm(x1);
fprintf('Diff. of last two: %6.2e\n\n',dx12)

%
figure; set(gca,'fontsize',12)
plot(1:n,x1,'b.',1:n,x2,'ro',1:n,xs,'k*');
legend('YALL1','SPGL1','Exact')
%}