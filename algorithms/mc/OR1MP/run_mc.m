function M_hat = run_mc(params)
  M = params.M;
  Idx = params.Idx;
  [numr,numc] = size(M);
  rank = 2;
  Known = find(Idx);
  data = M(Known);
  [U,Theta,V] = OR1MP(numr,numc,rank,Known,data);
  M_hat = U*diag(Theta)*V';
end
