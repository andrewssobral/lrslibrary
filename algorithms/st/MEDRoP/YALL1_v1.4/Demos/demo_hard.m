% Two ``hard'' problems
% prob 1: 150 nonzeros
% prob 2:   8 nonzeros

if ~exist('hard150.mat','file') || ~exist('hard8nz.mat','file')
    error('Data files hard150.mat and hard8nz.mat are not found. Please download them from the YALL1 website.');
end

% edit the following line to access SPGL1
addpath(genpath('/Users/yin/Documents/matlab/spgl1-1.7/'));
% edit the following line to access hard150 and hard8nz
addpath(genpath('/Users/yin/1st/matlab/YALL1/working/TestProbs/'));

clear
for prob = 1:2
    
switch prob
    case 1;
        load hard150;
    case 2; 
        load hard8nz;
        [Q,R] = qr(A',0);
        A = Q'; b = R'\b;
        clear Q R;
end

[m,n] = size(A); k = sum(xs ~= 0);
fprintf('\n*** Prob %i: size [n,m,k] = [%i,%i,%i]\n',prob,n,m,k);

% call yall1 solver
disp('--- YALL1 ---');
opts.tol = 1e-5;
opts.rho = 1e-5; 
opts.maxit = 3000;
opts.mu = 1/prob;
tic; [x,Out] = yall1(A, b, opts); toc;
rerr = norm(x-xs)/norm(xs);
fprintf('Iter %i: Rel_err = %6.2e\n',Out.iter,rerr);
fprintf('||x||_1 = %g\n',norm(x,1));
%close all; plot(1:n,xs,'ro',1:n,x,'b.'); drawnow

if exist('spgl1','file')
    disp('--- SPGL1 ---');
    spg_opts = spgSetParms('verbosity',0,'bpTol',opts.tol);
    tic; [x,r,g,info] = spgl1(A, b, 0, 0, [], spg_opts); toc
    rerr = norm(x-xs)/norm(xs);
    fprintf('[nA,nAt]=[%i,%i]:  Rel_err = %6.2e\n',...
        info.nProdA,info.nProdAt,rerr);
    fprintf('||x||_1 = %g\n',norm(x,1));
end
fprintf('\n')

end
