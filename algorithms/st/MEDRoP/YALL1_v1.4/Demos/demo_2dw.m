function demo_2dw(imgsize)

% Tests on 2D phantom images with wavelet basis
% Requires Image Processing  and wavelet Toolboxes

if nargin == 0; imgsize = 256; end

xsize = [imgsize imgsize];
n = prod(xsize);
m = floor(.5*n);
sigma = 0.005;

% generate xs
ph = phantom(imgsize);
xs = ph(:);
k = sum(xs ~= 0);
fprintf('\nSize [n,m,k] = [%i,%i,%i]\n',n,m,k);

p = randperm(n);
picks = sort(p(1:m)); picks(1) = 1;
noise = sigma*randn(m,1);
perm = randperm(n);

% set options
digit = 4; if sigma > 0; digit = 2; end 
opts.tol = 2*10^(-digit);

if exist('wavedec2','file')
    opts.basis = wavelet_basis(xsize);
    Mtype = {' pdct','pdwht'};
else
    disp('Wavelet not available, using DCT.')
    opts.basis.times = @(x)reshape( dct2(reshape(x,xsize)),numel(xs),1);
    opts.basis.trans = @(x)reshape(idct2(reshape(x,xsize)),numel(xs),1);
    Mtype = {' pdft','pdwht'};
end

for j = 1:2
    
    opts.rho = (j-1)*1e-3;
    fprintf('--- rho = %g ---\n',opts.rho)
    figure(j)
    
    % Original Image
    
    subplot(1,3,1); set(gca,'fontsize',16);
    imshow(reshape(xs,xsize),[]); title('Original');
    xlabel([int2str(imgsize) ' x ' int2str(imgsize)]);
    
    t0 = cputime;
    for k = 1:length(Mtype)
        op_A = eval(['@' Mtype{k} '_operator']);
        A = feval(op_A,picks,perm);
        b = A.times(xs) + noise;
        [x,Out] = yall1(A,b,opts);
        x = real(x); rerr = norm(x-xs)/norm(xs);
        Mat = upper(Mtype{k});
        fprintf([Mat ': iter %i  Rel_err = %6.2e\n'],...
            Out.iter,rerr);
        subplot(1,3,1+k); set(gca,'fontsize',16);
        imshow(reshape(x,xsize),[]); title(Mat);
        xlabel([int2str(m/n*100) '% data']);
    end
    fprintf('*CPU*: %6.2f sec.\n',cputime-t0);

end
