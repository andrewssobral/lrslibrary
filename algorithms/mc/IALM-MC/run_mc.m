function M_hat = run_mc(params)
  M = params.M;
  A = inexact_alm_mc(M,1e-4,100);
  M_hat = A.U*A.V';
end
