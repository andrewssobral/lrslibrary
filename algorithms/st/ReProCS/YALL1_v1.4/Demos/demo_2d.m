function demo_2d(imgsize)

% Tests on 2D data: phantom images
% Requires Image Processing Toolbox

if nargin == 0; imgsize = 256; end

xsize = [imgsize imgsize];
n = prod(xsize);
m = floor(.66*n);

% generate xs
Img = phantom(imgsize);
xs = Img(:);
k = sum(xs ~= 0);
fprintf('\nSize [n,m,k] = [%i,%i,%i]\n',n,m,k);

% noise and picks
sigma = 0.00;
noise = sigma*randn(m,1);
p = randperm(n);
picks = sort(p(1:m)); picks(1) = 1;
perm = randperm(n);

% set options
digit = 3; 
opts.tol = 2*10^(-digit);
Mtype = {' pdct','pdwht',' pdft'};

%opts.nonorth = 1;
for j = 1:2
    
    figure(j)
    opts.rho = (j-1)*1e-5;
    subplot(221); set(gca,'fontsize',16);
    imshow(reshape(xs,xsize),[]); title('Original');
    
    t0 = cputime;
    for k = 1:3
        op_handle = ['@' Mtype{k} '_operator'];
        A = feval(eval(op_handle),picks,perm);
        b = A.times(xs) + noise;
        [x,Out] = yall1(A, b, opts);
        x = real(x); rerr = norm(x-xs)/norm(xs);
        Mat = upper(Mtype{k});
        fprintf([Mat ': iter %3i  Rel_err = %6.2e\n'],...
            Out.iter,rerr);
        subplot(221+k); set(gca,'fontsize',16);
        imshow(reshape(x,xsize),[]); title(Mat);
    end
    fprintf('  rho = %4.1e,  CPU %6.2f sec.\n\n',...
        opts.rho,cputime-t0);
    
end
