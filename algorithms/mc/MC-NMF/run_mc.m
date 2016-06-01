function M_hat = run_mc(params)
  M = params.M;
  Idx = params.Idx;
  
  Omega = find(Idx);
  A = M(Omega);
  
  [m,n] = size(M);
  r = 2; t = 1;
  esr = ceil(t*r); 

  opts.tol = 1e-5;
  opts.maxit = 500;
  opts.print = 1;

  [X,Y,Out] = mc_nmf(A,Omega,esr,m,n,opts);
  M_hat = X*Y;
end
