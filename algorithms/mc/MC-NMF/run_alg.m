[numr,numc] = size(M);
Idx = randi([0 1],numr,numc); % ones(size(M));

params.M = M;
params.Idx = Idx;

L = run_mc(params);
S = (M - L);
