[numr,numc] = size(M);
I = randi([0 1],numr,numc); % ones(size(M));

M(M == 0) = 1e-3;
params.M = M.*I;
params.Idx = I;

L = run_mc(params);
S = (M - L);
