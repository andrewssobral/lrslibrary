function f = cost(X)
  global P;
  global PA;
  % Note that it is very much inefficient to explicitly construct the
  % matrix X in this way. Seen as we only need to know the entries
  % of Xmat corresponding to the mask P, it would be far more
  % efficient to compute those only.
  Xmat = X.U*X.S*X.V';
  f = .5*norm( P.*Xmat - PA , 'fro')^2;
end