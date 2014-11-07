% test l1 filtering

clear;
clc

% generating data
m = 1000;
n = 1000;
rho_s = 0.01;
rho_r = 0.01;
r = round(rho_r*min(m,n)) ;
U = (randn(m,r)); V = (randn(n,r));
A0 = U*V' ;

p = rho_r;
temp = randperm(m*n) ;
numCorruptedEntries = round(p*m*n) ;
corruptedPositions = temp(1:numCorruptedEntries) ;
E0 = zeros(m,n) ;
E0(corruptedPositions) = 100 * (rand(numCorruptedEntries,1)-0.5) ;

D = A0 + E0;

% ground truth

A0_rank = rank(A0,1e-3*norm(A0,2));
A0_nuclear = sum(svd(A0));
E0_0 = length(find(E0~=0));
E0_1 = sum((abs(E0(:))));

disp(['rank(A0)=', num2str(A0_rank), ', |A0|_*=', num2str(A0_nuclear),...
    ', |E0|_0=', num2str(E0_0), ', |E|_1=', num2str(E0_1)]);

% l1 filtering
[A_full, E_full, t_seed, t_l1f] = rpca_l1f(D);
% M = I0(1:1189,:); M = I0;
% clc; [A_full, E_full, t_seed, t_l1f] = rpca_l1f(M);
% imshow(E_full,[]);
t = t_seed + t_l1f;
disp(['Total Time:', num2str(t), ', Seed Time:',...
    num2str(t_seed), ', Filtering Time:', num2str(t_l1f)]);

A_rank = rank(A_full,1e-3*norm(A_full,2));
A_nuclear = sum(svd(A_full));
disp(['rank(A)=', num2str(A_rank), ', |A|_*=', num2str(A_nuclear)]);

eps = 1e-5;
idx = find(abs(E_full)<eps);
E_full(idx) = 0;
E_0 = length(find(E_full~=0));
E_1 = sum((abs(E_full(:))));
disp(['|E|_0=', num2str(E_0), ', |E|_1=', num2str(E_1)]);

errA = norm(A_full-A0,'fro')/norm(A0,'fro');
dif = A0-A_full;
DIF = max(abs(dif(:)));
AVDIF = mean(abs(dif(:)));
disp(['error:', num2str(errA)]);

