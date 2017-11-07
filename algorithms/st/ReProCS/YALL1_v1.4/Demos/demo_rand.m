function demo_rand(n)

% Test on rand problems with Gaussian matrices
% Data can contain white or impulsive noise or both

% uncomment and edit the following line to run SPGL1
%addpath(genpath('/Users/yin/Documents/matlab/spgl1-1.7/'));

if nargin == 0; n = 1000; end
m = floor(0.3*n);
k = floor(0.2*m);
fprintf('\nSize [n,m,k] = [%i,%i,%i]\n',n,m,k);

% generate A, b, xs
sigma = 0.00;  % white noise std
perc = 0;      % impulsive noise percentage
nonorth = 1;
[A,b,xs] = gen_data(m,n,k,sigma,perc,nonorth);

% call solver
digit = 6; if sigma > 0; digit = 3; end
opts.tol = 5*10^(-digit);
opts.print = 0;
opts.nonorth = nonorth;
disp('--- YALL1:  BP');
opts.nu = 0;
opts.rho = 0;
tic; [x,Out] = yall1(A,b,opts); toc
rerr = norm(x-xs)/norm(xs);
fprintf('Iter %4i: Rel_err = %6.2e\n\n',Out.iter,rerr);

disp('--- YALL1: L1/L1');
opts.nu = .5;
opts.rho = 0;
tic; [x,Out] = yall1(A,b,opts); toc
rerr = norm(x-xs)/norm(xs);
fprintf('Iter %4i: Rel_err = %6.2e\n\n',Out.iter,rerr);

disp('--- YALL1: L1/L2');
opts.nu = 0;
opts.rho = 1e-3; 
tic; [x,Out] = yall1(A,b,opts); toc
rerr = norm(x-xs)/norm(xs);
fprintf('Iter %4i: Rel_err = %6.2e\n\n',Out.iter,rerr);

if exist('spgl1','file')
    disp('--- SPGL1:');
    spg_opts = spgSetParms('verbosity',0,'bpTol',opts.tol);
    tic; [x,~,~,info] = spgl1(A, b, 0, 0, [], spg_opts); toc
    rerr = norm(x-xs)/norm(xs);
    fprintf('[nA,nAt]=[%i,%i]:  Rel_err = %6.2e\n\n',...
        info.nProdA,info.nProdAt,rerr);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,b,xs] = gen_data(m,n,k,sigma,perc,nonorth)

A = randn(m,n);
d = 1./sqrt(sum(A.^2));
A = A*sparse(1:n,1:n,d);
xs = zeros(n,1);
p = randperm(n);
xs(p(1:k)) = randn(k,1);
b = A*xs;
if ~nonorth
    [Q,R] = qr(A',0);
    A = Q'; b = R'\b; 
end

% white noise
white = sigma*randn(m,1);
b = b + white;

% impulsive noise
p = randperm(m);
L = floor(m*perc/100);
b(p(1:L)) = 2*(rand(L,1) > 0) - 1;
