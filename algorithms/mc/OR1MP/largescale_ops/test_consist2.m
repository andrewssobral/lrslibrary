function test_consist2()
  S = sparse(diag(eye(5)));
  idx = find(S~=0);

  v = rand(5, 1);

  out_us_nat = sparse_update_native(S, v);
  out_us_mex = S; 
  sparse_update(out_us_mex, v);

  if nnz(out_us_mex ~= out_us_nat)>0
      warning('inconsistent results found in sparse_update')
  else
      disp('[MEX: sparse_update] consistency check passed.');
  end
end