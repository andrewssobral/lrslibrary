% --------------------------------------------
% demo_rand: random low-rank and sparse
% --------------------------------------------
% This demo tests the matlab code lmafit_sms.m
% for solving sparse matrix separation problem.
% July 25, 2010. Yin Zhang

function demo_rand(Solver,m,n,r,impulse,tol)

if nargin < 1; Solver = 1; end
if nargin < 2; m = 400; end
if nargin < 3; n = m; end
if nargin < 4; r = round(.1*min(m,n)); end
if nargin < 5; impulse = 0.1; end 
if nargin < 6; tol = 1e-4; end
    
if ~exist('PROPACK','dir');
   path(path,genpath(pwd));
   warning off all;
end

sigma = 0e-3; % noise STD
tol = max(tol, .1*sigma); 
tol = min(1e-1,tol); 
maxit = 200;

% lmafit parameter setup
opts.tol = tol;
opts.maxit = maxit;
opts.est_rank = 1;
opts.print = 0;
% initial rank estimate
k = round(0.25*min(m,n));

% random low-rank and sparse
XY0 = randn(m,r)*randn(r,n);            % random low-rank
Z0 = rand(m,n) < impulse;               % random sparse pattern
Z0 = Z0.*randn(m,n);                    % Gaussian values
%disp(max(abs(XY0(:)))/max(abs(Z0(:))))

% data matrix D (with possible noise)
D = XY0 + Z0;
D = D + sigma*randn(m,n);

beta = []; %beta = 1.2.^(0:20);
t = 1;

fprintf('\n[m n r k] = [%i, %i, %i, %i]\n',m,n,r,k);
fprintf('noise sigma: %6.2e\n',sigma);
fprintf('---------------------------------\n')

% call solvers
for sol = Solver
    tic
    switch sol
        case 1; fprintf('   ... LMafit_sms ...\n')
            [X,Y,Z,out] = lmafit_sms_v1(D,k,opts,beta);
            XY = X*Y; iter = out.iter; rk = size(X,2);
        case 2; fprintf('   ... IALM_rpca  ...\n')
            lambda = t/sqrt(m);
            fprintf('lambda = %.2f/sqrt(m)\n',t);
            [XY,Z,iter] = inexact_alm_rpca(D,lambda,tol/10,maxit);
            rk = rank(XY);
    end
    toc
    % error calculation function
    ferr = @(U,V) norm(U-V,'fro')/norm(V,'fro');
    fprintf('XY RelErr = %8.3e\n',ferr(XY,XY0));
    fprintf('Z  RelErr = %8.3e\n',ferr( Z,Z0 ));
    fprintf('A  RelErr = %8.3e\n',ferr(XY+Z,D));
    fprintf('iter %i: final rank = %i\n',iter,rk);
    fprintf('---------------------------------\n')
end
fprintf('\n\n')