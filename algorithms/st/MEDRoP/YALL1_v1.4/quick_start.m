clear;
if ~exist('Utilities','dir'); Run_Me_1st; end

% problem sizes 
n = 1000; m = 300; k = 60; 
sigma = 0.01;
opts.rho = eps;

nrun = 10;
Iter = zeros(nrun,1);
Err = zeros(nrun,1);
T = zeros(nrun,1);

for j = 1:nrun
    
    % generate (A,xs,b)
    A = randn(m,n);
    xs = zeros(n,1);
    p = randperm(n);
    xs(p(1:k)) = randn(k,1);
    b = A*xs + sigma*randn(m,1);
    
    % (orth)normalize the rows of A
    if ~exist('nonorth','var');
        nonorth = randn > 0;
    end
    
    if nonorth;
        d = 1./sqrt(sum(A.^2,2));
        A = sparse(1:m,1:m,d)*A;
        b = d.*b;
    else
        [Q, R] = qr(A',0);
        A = Q'; b = R'\b;
    end
    
    % call YALL1
    opts.tol = 5e-8;
    if sigma > 0;
        opts.tol = 5e-3;
        opts.rho = sigma;
    end
    opts.print = 0;
    t0 = tic; [x,Out] = yall1(A, b, opts);
    err = norm(x-xs)/norm(xs);
    fprintf('nonorth: %i, iter = %4i, relative error = %e\n',...
        nonorth,Out.iter,err)
    Iter(j) = Out.iter;
    Err(j) = err;
    T(j) = toc(t0);
end

fprintf('\n[n,m,k] = [%i,%i,%i], %i runs\n',n,m,k,nrun);
fprintf('Average: iter %3i, rel_err %6.2e, time %6.2e\n\n',...
    round(mean(Iter)),mean(Err),mean(T))
