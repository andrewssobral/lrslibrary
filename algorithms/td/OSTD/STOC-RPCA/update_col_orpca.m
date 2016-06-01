function L = update_col_orpca(L,A,B,lambda1)
  [junk,ncol] = size(L);
  A = A + lambda1*eye(ncol,ncol);
  for j = 1:ncol
    bj = B(:,j);
    lj = L(:,j);
    aj = A(:,j);
    temp = (bj-L*aj)/A(j,j) + lj;
    L(:,j) = temp/max(norm(temp),1);
  end
end