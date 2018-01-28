clear

n = 64^2; % signal length
m = n/4;  % number of measurements
k = m/8;  % number of nonzeros

Dyna = 100:50:300;

for dyna = Dyna
    
    % random signal with dyna-dB dynamic range
    valx = dyna/20*rand(k,1);
    valx = valx - min(valx);
    valx = valx/max(valx)*dyna/20;
    xs = zeros(n,1); p = randperm(n);
    xs(p(1:k)) = 10.^valx.*sign(randn(k,1));
    
    l1_opt = norm(xs,1);
    fprintf('dyna = %g, |x*|_1 = %6.3e\n',dyna,l1_opt);
    
    % partial dct measurement matrix
    Omega = randperm(n);
    Omega = Omega(1:m);
    Omega = sort(Omega);
    op_A = @pdct_operator;
    A = feval(op_A,Omega,1:n);
    
    % observations
    b = A.times(xs);
    
    opts.tol = eps; opts.maxit = 1000;    
    tic, [x,Out] = yall1(A, b, opts); toc
    fprintf('Iter: %i\n',Out.iter)
    fprintf('Abs. Ax-b  Residual: %6.3e\n',norm(A.times(x)-b));
    fprintf('Abs. solution error: %6.3e\n',norm(x-xs));
    fprintf('Rel. Ax-b  Residual: %6.3e\n',norm(A.times(x)-b)/norm(b));
    fprintf('Rel. solution error: %6.3e\n\n',norm(x-xs)/norm(xs));
    
end
