function M_hat = run_mc(params)

  avg = mean(mean(params.M, 1), 2);
  M = params.M - avg;

  [U_t, SV_t] = ncrmc(M, params.Idx);
  M_hat = U_t * SV_t + avg;

end
