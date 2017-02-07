function [ X ] = X_Solver_first(D,rho)
  [U,S,V] = svd(D);
  S0 = diag(S);
  r = length(S0);
  P = [ones(r,1), 1-S0, 1/2/rho-S0];
  rt = zeros(r,1);
  for t = 1:r
    p = P(t,:);
    Delta = p(2)^2-4*p(1)*p(3);
    if Delta <= 0
      rt(t) = 0;
    else
      rts = roots(p);
      rts = sort(rts);
      if rts(1)*rts(2)<=0
        rt(t) = rts(2);
      elseif rts(2)<0
        rt(t) = 0;
      else
        funval = log(1+rts(2))+rho*(rts(2)-S0(t)).^2;
        if funval>log(1+0)+rho*(0-S0(t)).^2;
          rt(t) = 0;
        end
      end
    end
  end

  SSS = diag(rt);
  [m,n] = size(D);
  sig = zeros(m,n);
  sig(1:min(m,n),1:min(m,n)) = SSS;

  X = U*sig*V';
end