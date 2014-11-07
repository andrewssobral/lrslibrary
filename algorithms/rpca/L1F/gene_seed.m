
function [A_seed, column_seed, row_seed, t_seed] = gene_seed(D_full, sr, sc)

t_seed = 0;

if nargin < 3
    sc = 10;
end
if nargin < 2
    sr = 10;
end

[m, n] = size(D_full);

d_r = ceil(m/sr);
d_c = ceil(n/sc);

row = randperm(m);
row_seed = row(1:d_r);
column = randperm(n);
column_seed = column(1:d_c);

D_seed = D_full(row_seed, column_seed);
%[A_seed, ~] = ialm_rpca(D_seed);
[A_seed, ~] = inexact_alm_rpca(D_seed);
rA = rank(A_seed);

%% estimate the rank
while d_r/rA < sr || d_c/rA < sc || max(m/d_r, n/d_c) < 0.5
        d_r = sr*rA;
        d_c = sc*rA;
        
        row = randperm(m);
        clc;
        row
        d_r
        row_seed = row(1:d_r);
        column = randperm(n);
        column_seed = column(1:d_c);
        
        D_seed = D_full(row_seed, column_seed);
        %[A_seed, ~] = ialm_rpca(D_seed);
        [A_seed, ~] = inexact_alm_rpca(D_seed);
        rA = rank(A_seed);
end
true_r = rA;

% prior strategy
% S = svd(A_seed);
% ds = diag(S);
% if rA/min(d_r, d_c) > 0.25
%     sr = ds(1:rA-1)./ds(2:rA);
%     [~, true_r] = max(sr(1:end-1));
% else
%     true_r = rA;
% end

%% seed recovery
d_r = sr*true_r;
d_c = sr*true_r;
row = randperm(m);
row_seed = row(1:d_r);
column = randperm(n);
column_seed = column(1:d_c);

D_seed = D_full(row_seed, column_seed);

%ts = tic;
%[A_seed, ~] = ialm_rpca(D_seed);
[A_seed, ~] = inexact_alm_rpca(D_seed);
%t_seed = toc(ts);

