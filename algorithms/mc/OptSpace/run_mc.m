function M_hat = run_mc(params)
  M = params.M;
  tol = 1e-8;
  [X,S,Y,dist] = OptSpace(sparse(M),[],20,tol);
  M_hat = X*S*Y';
end
