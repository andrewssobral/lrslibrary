function S = Tproj_partial(S, a_col, a_row)
% The T() operator for sparse prejection

[d1, d2] = size(S);
kcol = floor(a_col*d1);
krow = floor(a_row*d2);
[~, colloc] = maxk(abs(S), kcol, 1, 'sorting', false);
[~, rowloc] = maxk(abs(S'), krow, 1, 'sorting', false);

%[iS, jS, vS] = find(S);

Jid = repmat((1:d2), kcol, 1); Jid = Jid(:);
Iid = repmat((1:d1), krow, 1); Iid = Iid(:);

%rowlocT = rowloc';
mask_col = sparse(colloc(:), Jid, 1, d1, d2, nnz(S));
mask_row = sparse(Iid, rowloc(:), 1, d1, d2, nnz(S));

S = S.*mask_col.*mask_row;

