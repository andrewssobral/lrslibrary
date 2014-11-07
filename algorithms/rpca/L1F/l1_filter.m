function [A_full, E_full, t_l1f] = l1_filter(D_full, column_seed, row_seed, A_seed)

t_l1f = 0;

[m, n] = size(D_full);
if nargin < 4
    D_seed = D_full(row_seed, column_seed);
    %[A_seed, ~] = ialm_rpca(D_seed);
    [A_seed, ~] = inexact_alm_rpca(D_seed);
end

rA = rank(A_seed);

[U,S,V] = svd(A_seed);
Ur = U(:,1:rA);
Vr = V(:,1:rA);
Sr = diag(S);
Sr = diag(Sr(1:rA));

%ts = tic;

As_row = Ur*Sr;
As_column = Sr*Vr';

A_row = zeros(rA, n);
A_column = zeros(rA, m);

column_comp = setdiff(1:n, column_seed);
pinvUr = Ur';
[~, Ac_row] = solve_ml1(D_full(row_seed, column_comp), Ur, pinvUr);
A_row(:, column_comp) = Ac_row;
A_row(:, column_seed) = As_column;

row_comp = setdiff(1:m, row_seed);
pinvVr = Vr';
[~, Ac_column] = solve_ml1(D_full(row_comp, column_seed)', Vr, pinvVr);
A_column(:, row_comp) = Ac_column;
A_column(:, row_seed) = As_row';

A_full = A_column'*inv(Sr)*A_row;
E_full = D_full - A_full;

%t_l1f = toc(ts);