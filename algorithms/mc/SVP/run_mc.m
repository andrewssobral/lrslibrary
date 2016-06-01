function M_hat = run_mc(params)
  M = params.M;
  Idx = params.Idx;
  warning('off','all');
  [numr,numc] = size(M);
  rank = 2;
  Known = find(Idx);
  data = M(Known);
  [U,S,V] = svp(Known,data,numr,numc,rank);
  M_hat = U*diag(S)*V';
  warning('on','all');
end
