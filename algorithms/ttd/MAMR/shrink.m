function Y = shrink(X,a)
% Copyright:Xinchen YE, Tianjin University, 2014
  tmpY = max(X-a,0);
  Y = tmpY + min(X+a,0);
end