function G = egrad(X)
  global P;
  global PA;
  % Same comment here about Xmat.
  Xmat = X.U*X.S*X.V';
  G = P.*Xmat - PA;
end