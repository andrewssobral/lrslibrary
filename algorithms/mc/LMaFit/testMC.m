% script for tesing the matrix completion code lmafit_mc

clear

% problem parameters
% disp('Example Input: [5000 5000 10 6]') 
% params = input('[m n r dr] = ');
% 
% m  = params(1);
% n  = params(2);
% r  = params(3); 
% dr = params(4);

m = 500;  % row of the matrix
n = 500;  % column of the matrix
r = 5;    % rank of the matrix
dr = 6;    % degree of freedom of the matrix

if dr <= 1; error('dr > 1 is required'); end

% sample ratio: sr
sr = min(dr*r*(m+n-r)/(m*n),1);
if sr >= 1; error('sample ratio sr >= 1'); end

% rank estimate
k = ceil(1.25*r);

% noise
sigma = 0e-2;
fprintf('k = %i, noise sigma: %6.2e\n',k,sigma);

% generate problem data
fprintf('Generating ...')
RandStream.setDefaultStream...
     (RandStream('mt19937ar','seed',sum(100*clock)));
Xs = randn(m,r);
Ys = randn(r,n);
L = ceil(sr*m*n);

Known = unique(round(rand(1,2*L)*m*n));
if Known(1) == 0; Known = Known(2:end); end
Lc = length(Known);
if Lc > L
    p = randperm(Lc);
    Known = sort(Known(p(1:L)));
else
    L = Lc;
end

% generate data
[Ik, Jk] = ind2sub([m n],Known);
data = partXY(Xs',Ys,Ik,Jk,L);

% add noise
noise = randn(size(data));
level = sigma*norm(data)/norm(noise);
data = data + level*noise; 
fprintf(' Done. sr = %6.4f\n',L/(m*n));
clear noise;

% error calculation functions

lrdot = @(X1,Y1,X2,Y2) sum(sum((X1'*X2).*(Y1*Y2')));

% problem specification
opts = [];
%opts.tol = 2e-4;
%opts.maxit = 200;
%opts.Zfull = 0;
%opts.DoQR  = 1;
% opts.print = 0;
%opts.est_rank = 1;

% call solver
tic; [X,Y,Out] = lmafit_mc_adp(m,n,k,Known,data,opts); toc
clear Prob;

% compute relerr for the data part
comp = partXY(X',Y,Ik,Jk,L);
ferr = @(u,v) norm(u(:)-v(:))/norm(v(:));
fprintf('RelErr_Known = %8.3e\n',ferr(comp, data));


% compute full relerr = ||XsYs - XY||_F/||XsYs||_F 
XsYsnrm2 = lrdot(Xs,Ys,Xs,Ys);   % ||XsYs||_F^2
XYnrm2 = lrdot(X,Y,X,Y);         % ||XY||_F^2
cross = lrdot(Xs,Ys,X,Y);        % <XsYs, XY>
rerr2 = (XsYsnrm2 + XYnrm2 - 2*cross)/XsYsnrm2;
fprintf('RelErr_Full  = %8.3e\n',sqrt(rerr2));
fprintf('Est_rank  = %i\n\n',Out.rank);

if 0
    semilogy(1:Out.iter,Out.obj,'b.-')
    set(gca,'fontsize',16);
    xlabel('Iteration'); ylabel('Objective');
end
