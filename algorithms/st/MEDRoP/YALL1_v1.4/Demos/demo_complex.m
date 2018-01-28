function demo_complex(n)

% test solver for complex BP

close all;

opts.nonorth = 1;
if nargin == 0; 
    if opts.nonorth; n = 1200; else n = 1024*16; end
end
m = floor(n/4); k = floor(m/4);
fprintf('\nSize [n,m,k] = [%i,%i,%i]\n',n,m,k);

% generate xs
xs = zeros(n,1);
p = randperm(n); p = p(1:k);
xs(p) = randn(k,1) + 1i*randn(k,1);

% generate partial DFT data
p = randperm(n);
picks = sort(p(1:m)); picks(1) = 1;
perm = 1:n; % column permutations allowable

if opts.nonorth
    A = randn(m,n) + 1i*randn(m,n);
    d = 1./sqrt(sum(abs(A).^2));
    A = A*sparse(1:n,1:n,d);
    b = A*xs;
else
    A = pdft_operator(picks,perm);
    b = A.times(xs);
end

% add noise
sigma = norm(b,inf)*0.001;  % noise std
noise = randn(m,1) + 1i*randn(m,1);
b = b + sigma*noise;

% set solver options
digit = 6; if sigma > 0; digit = 3; end
opts.tol = 5*10^(-digit);
plotting = 1;

for j = 1:2
    
    opts.rho = (j-1)*1e-3;
    fprintf('--- rho = %g ---\n',opts.rho);
    % call solver
    tic; [x,Out] = yall1(A, b, opts); toc
    rerr = norm(x-xs)/norm(xs);
    fprintf('iter %4i  Rel_err = %6.2e\n',Out.iter,rerr);
    
    % plot
    if plotting
        figure(j);
        subplot(211); set(gca,'fontsize',16)
        plot(1:n,real(xs),'ro',1:n,real(x),'b.');
        title('Real Part')
        subplot(212); set(gca,'fontsize',16)
        plot(1:n,imag(xs),'ro',1:n,imag(x),'b.');
        title('Imag Part')
    end
    
end
