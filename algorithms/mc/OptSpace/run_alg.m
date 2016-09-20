% MC | OptSpace | A Matrix Completion Algorithm (Keshavan et al. 2009)
% process_video('MC', 'OptSpace', 'dataset/demo.avi', 'output/demo_MC-OptSpace.avi');

numr = size(M,1);
numc = size(M,2);
Idx = randi([0 1],numr,numc);
MI = M.*Idx; % show_2dvideo(MI,m,n);
tol = 1e-8;
[X,S,Y,dist] = OptSpace(sparse(MI),[],20,tol);
L = X*S*Y'; % show_2dvideo(L,m,n);
S = M - L;
