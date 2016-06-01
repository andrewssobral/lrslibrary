function test_consist1()
  U = rand(500, 3);
  V = rand(3, 10000); 

  iR = [1, 3, 3];
  iC = [2, 4, 6];

  len = length(iR);


  out_us_mex = sparse_inp(U', V, iR, iC)';
  out_us_net = sparse_inp_native(U, V, iR, iC); % call native 

  if nnz(out_us_mex ~= out_us_net)>0
      warning('inconsistent results found in sparse_inp_native')
  else
      disp('[MEX: sparse_inp] consistency check passed.');
  end
end