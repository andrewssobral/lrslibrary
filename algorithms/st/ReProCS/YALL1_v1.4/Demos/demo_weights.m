function demo_weights(n)

% Test weighted 1-norm on rand problems 
% with Gaussian matrices

if nargin == 0; n = 1000; end

m = floor(0.3*n);
k = floor(0.4*m);
fprintf('\nSize [n,m,k] = [%i,%i,%i]\n',n,m,k);

% generate A, b, xs
sigma = 0.001;  % noise std
opts.nonorth = 1;
[A,b,xs,w] = gen_data(m,n,k,sigma,opts.nonorth);

% set options
digit = 6; if sigma > 0; digit = 3; end
opts.tol = 10^(-digit);
opts.weights = w;
opts.print = 0;

% call solver
disp('--- BP:')
tic; [x,Out] = yall1(A,b,opts); toc
rerr = norm(x-xs)/norm(xs);
fprintf('Iter %4i: Rel_err = %6.2e\n',Out.iter,rerr);

disp('--- L1L1:')
opts.nu = .5;
tic; [x,Out] = yall1(A,b,opts); toc
rerr = norm(x-xs)/norm(xs);
fprintf('Iter %4i: Rel_err = %6.2e\n',Out.iter,rerr);

disp('--- L1L2:')
opts.nu = 0; opts.rho = 1e-3;
tic; [x,Out] = yall1(A,b,opts); toc
rerr = norm(x-xs)/norm(xs);
fprintf('Iter %4i: Rel_err = %6.2e\n',Out.iter,rerr);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,b,xs,w] = gen_data(m,n,k,sigma,nonorth)

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
b = b + sigma*randn(m,1);

% weights
w = ones(n,1);
pk = randperm(k);
iw0 = p(pk(1:floor(k/2)));
w(iw0) = 0; 