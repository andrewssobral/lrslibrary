% --------------------------------------------
% demo_rand: random low-rank and sparse
% --------------------------------------------
% This demo tests the matlab code lmafit_sms.m
% for solving sparse matrix separation problem.
% July 25, 2010. Yin Zhang

function demo_chkb(Solver,m,impulse)

if nargin < 1; Solver = 1; end
if nargin < 2; m = 256; end
if nargin < 3;
    disp('impulse [0,1]')
    impulse = input(' impulse = ');
end

if ~exist('PROPACK','dir');
    path(path,genpath(pwd));
    warning off all;
end

n = m; r = 2;
sigma = 0e-3; % noise STD
tol = max(2e-4, 0.5*sigma);
maxit = 200;

% lmafit parameter setup
opts.tol = tol;
%opts.maxit = maxit;
%opts.est_rank = true;
%opts.min_rank = 1;
beta = 1:10;

% random low-rank and sparse
XY0 = checkerboard(m/8);                % checkerboard
Z0 = rand(m,n) < impulse;               % random sparse pattern
Z0 = Z0.*randn(m,n);                    % Gaussian values

% data matrix D (with possible noise)
D = XY0 + Z0;
D = D + sigma*randn(m,n);
k = round(m/2);

subplot(1,1+numel(Solver),1); imshow(D,[]);
xlabel(sprintf('size: %i x %i',m,n),'fontsize',14);
title(sprintf('original'),'fontsize',14); drawnow;
 
% initial rank estimate
fprintf('\n[m n r k] = [%i, %i, %i, %i]\n',m,n,r,k);
fprintf('noise sigma: %6.2e\n',sigma);
fprintf('---------------------------------\n')

SolverName = ['LMaFit sms'; 'IALM rpca ']; 
% call solvers
for sol = Solver
    tic
    switch sol
        case 1; fprintf('   ... LMafit_sms ...\n')
            [X,Y,Z,out] = lmafit_sms_v1(D,k,opts,beta);
            XY = X*Y; iter = out.iter;
        case 2; fprintf('   ... IALM_rpca  ...\n')
            t = 0.5; lambda = t/sqrt(m);
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
    pos = 2+(sol-1)*(numel(Solver)-1);
    subplot(1,1+numel(Solver),pos); imshow(XY,[]);
    xlabel(sprintf('Rel-Err: %.2e',ferr(XY,XY0)),'fontsize',14);
    title(sprintf('%10s',SolverName(sol,:)),'fontsize',14); drawnow;

end
fprintf('\n\n')
