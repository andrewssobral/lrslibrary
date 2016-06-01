% TTD | 3WD | 3-Way-Decomposition (Oreifej et al. 2012)
% process_video('TTD', '3WD', 'dataset/demo.avi', 'output/demo_3WD.avi');
%{
clc;
lrs_load_conf;
load(fullfile(lrs_conf.lrs_dir,'dataset','trafficdb','traffic_patches.mat'));
V = im2double(imgdb{100});
[M,m,n,p] = convert_video3d_to_2d(V);
%}

D = M;

%CM = zeros(size(D));
[numr,numc] = size(D);
CM = randi([0 1],numr,numc); % simulated confidence map (binary matrix)

tauc = 100*.04; % .5 for with Pi outside / 5 for with Pi inside. lambda = lambdac/sqrt(m) 1.1
tau = tauc/sqrt(size(D,1));
lambdac = 1000; % 2000 good for Forb alone
lambda = lambdac/sqrt(size(D,1)) ;
tol = 1e-6;
maxIter = 100;
[A_dual,E_dual,O_dual,iter,Y] = ThreeWayDec(D,CM,tau,lambda,tol,maxIter);

% A(original) = L(low-rank) + S(sparse) + E(error)
L = A_dual;
S = O_dual;
E = E_dual;

M_hat = L + S + E;
error = norm(M_hat(:)-M(:))/norm(M(:));
disp(['Error: ' num2str(error)]);

%{
show_2dvideo(M,m,n);
show_2dvideo(A_dual,m,n);
show_2dvideo(E_dual,m,n);
show_2dvideo(O_dual,m,n);
%}