% --------------------------------------------
% demo_spim: SParse IMage (+ random low-rank)
% --------------------------------------------
% This demo tests the matlab code lmafit_sms.m
% for solving sparse matrix separation problem.
% July 25, 2010. Yin Zhang

function demo_spim(Solver,r,impulse)

if nargin < 1; Solver = 1:2; end
if nargin < 2; r = 10; end
if nargin < 3; impulse = 0.05; end

warning off all;
if ~exist('PROPACK','dir');
   path(path,genpath(pwd));
end

h = figure(9);
set(h,'units','normalized','outerposition',[0 .5 .8 .8]);

sigma = 0e-3; % noise STD
tol = max(2e-4,2*sigma);
maxit = 200;

% lmafit parameter setup
opts.tol = tol;
%opts.maxit = maxit;
%opts.est_rank = true;
%opts.min_rank = 1;

beta = 1:10;
Imgs = 1:5;

for ii = 1:numel(Imgs)

    imn = Imgs(ii);
    % input image
    switch imn
        case 1; I = zeros(329);
      I(29:300,:) = imread('blobs.png'); % make it square
        case 2; I = imread('circles.png');
        case 3; I = im2bw(phantom(256),.25);
        case 4; I = imread('text.png');
        case 5; I = imread('rice.png');
    end
    I = double(im2bw(I));

    % generate random low-rank
    [m,n] = size(I);
    XY0 = randn(m,r)*randn(r,n);            % random low-rank
    
    % impulsive noise + sparse pattern
    Z0 = rand(m,n) < impulse;               % (0,1) sparse
    %Z0 = Z0.*randn(m,n);                    % Gaussian sparse
    Z0 = min(Z0 + I,1);                     % add sparse pattern
    
    % data matrix D (with possible noise)
    D = XY0 + Z0;
    D = D + sigma*randn(m,n);

    % initial rank estimate
    k = min(5*r,round(0.25*min(m,n)));    
    fprintf('\n[m n r k] = [%i, %i, %i, %i]\n',m,n,r,k);
    fprintf('noise sigma: %6.2e\n',sigma);

    % call solvers
    for sol = Solver
        tic
        switch sol
            case 1; fprintf('   ... LMafit_sms ...\n')
                [X,Y,Z,out] = lmafit_sms_v1(D,k,opts,beta);
                XY = X*Y; iter = out.iter;
            case 2; fprintf('   ... IALM_rpca  ...\n')
                t = .5; lambda = t/sqrt(m);
                fprintf('lambda = %.2f/sqrt(m)\n',t);
                [XY,Z,iter] = inexact_alm_rpca(D,lambda,tol,maxit);
        end
        toc
        % error calculation function
        ferr = @(U,V) norm(U-V,'fro')/norm(V,'fro');
        fprintf('XY RelErr = %8.3e\n',ferr(XY,XY0));
        fprintf('Z  RelErr = %8.3e\n',ferr( Z,Z0 ));
        fprintf('A  RelErr = %8.3e\n',ferr(XY+Z,D));
        fprintf('iter %i: final rank = %i\n',iter,rank(XY));
        fprintf('---------------------------------\n')
        
        pos = (numel(Solver)-1)*(sol-1)*numel(Imgs)+ii;
        subplot(numel(Solver),numel(Imgs),pos); 
        imshow(Z,[]);
        xlabel(sprintf('Rel-Err: %.2e',ferr(Z,Z0)),...
            'fontsize',14);
        title(sprintf('Density: %.2f',nnz(Z)/(m*n)),...
            'fontsize',14);
        
    end
    drawnow;
    
end
fprintf('\n\n')