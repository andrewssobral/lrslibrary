[numr,numc] = size(M);
I = randi([0 1],numr,numc); % ones(size(M));

params.M = M;
params.Idx = I;

L = run_mc(params);
S = (M - L);
