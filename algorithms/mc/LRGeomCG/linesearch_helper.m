function t = linesearch_helper(X, H)
  % Note that you would not usually need the Hessian for this.
  residual_omega = nonzeros(problem.egrad(X));
  dir_omega      = nonzeros(problem.ehess(X, H));
  t = - dir_omega \ residual_omega ;
end