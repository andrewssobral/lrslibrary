function [U,S,V] = fastsvd(U,S,V,y,I,J,k)

  %Find the singular values of U S V' + Y
  %Y is a sparse matrix

  m = size(U,1);
  n = size(V,1);

  myAfunc = @(x) (U*(S*(V'*x)) + compOmegaYx(m,n,y,I,J,x));
  myAtfunc = @(x) (V*(S*(U'*x)) + compOmegaYtx(m,n,y,I,J,x));
  [U,S,V] = lansvd(myAfunc,myAtfunc, m,n,k,'L');

end
