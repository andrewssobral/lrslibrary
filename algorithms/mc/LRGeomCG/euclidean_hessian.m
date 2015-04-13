function ehess = euclidean_hessian(X, H)
  global P;
  global problem;
  % The function tangent2ambient transforms H (a tangent vector) into
  % its equivalent ambient vector representation. The output is a
  % structure with fields U, S, V such that U*S*V' is an mxn matrix
  % corresponding to the tangent vector H. Note that there are no
  % additional guarantees about U, S and V. In particular, U and V
  % are not orthonormal.
  ambient_H = problem.M.tangent2ambient(X, H);
  Xdot = ambient_H.U*ambient_H.S*ambient_H.V';
  % Same comment here about explicitly constructing the ambient
  % vector as an mxn matrix Xdot: we only need its entries
  % corresponding to the mask P, and this could be computed
  % efficiently.
  ehess = P.*Xdot;
end