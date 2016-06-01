function [ X,T ] = DC(D,rho,T0,a)
  [U,S,V] = svd(D,'econ');
  for t = 1:100  
    [ X,T1 ] = DCInner(S,rho,T0,a,U,V);
    err = sum((T1-T0).^2);
    if err < 1e-6
      break
    end
    T0 = T1;
  end
  T = T1;
end

function [ X,t ] = DCInner(S,rho,J,epislon,U,V)
  lambda=1/2/rho;
  S0 = diag(S);
  grad=(1+epislon)*epislon./(epislon+J).^2;
  t=max(S0-lambda*grad,0);
  X=U*diag(t)*V';
end
