[numr,numc] = size(M);
Idx = randi([0 1],numr,numc); % ones(size(M));

M(M == 0) = 1e-3;
params.M = M.*Idx;
params.Idx = Idx;

L = run_mc(params);
S = (M - L);
