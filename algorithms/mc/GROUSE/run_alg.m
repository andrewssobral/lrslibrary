% MC | GROUSE | Grassmannian Rank-One Update Subspace Estimation (Balzano et al. 2010)
% process_video('MC', 'GROUSE', 'dataset/demo.avi', 'output/demo_MC-GROUSE.avi');

[numr,numc] = size(M);
Idx = randi([0 1],numr,numc); % ones(size(M));
maxrank = 1;
maxCycles = 100;
step_size = 0.1;

[Usg, Vsg, err_reg] = grouse(M,Idx,numr,numc,maxrank,step_size,maxCycles);
L = Usg*Vsg';
S = M - L;

% show_2dvideo(M,m,n);
% show_2dvideo(M.*Idx,m,n);
% show_2dvideo(L,m,n);
% show_2dvideo(S,m,n);
