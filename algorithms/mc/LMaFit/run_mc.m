function M_hat = run_mc(params)
  M = params.M;
  Idx = params.Idx;
  [numr,numc] = size(M);
  rank = 10;
  Known = find(Idx);
  data = M(Known);
  [X,Y] = lmafit_mc_adp(numr,numc,rank,Known,data,[]);
  M_hat = X*Y;
end
