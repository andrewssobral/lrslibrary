% solve the problem:
%   min_{x,e} 0.5*|z-Dx-e|_2^2 + 0.5*lambda1*|x|_2^2 + lambda2*|e|_1
%
% solve the projection by APG
% input:
%   z - data point
%   D - basis matrix
%   lambda1, lambda2 - tradeoff parameters
% output:
%   r - projection coefficient
%   e - sparse noise
% copyright Jiashi Feng (jshfeng@gmail.com)
%
function [x,e] = solve_proj2(z,D,lambda1,lambda2)
  % initialization
  [ndim,ncol] = size(D);
  e = zeros(ndim,1);
  x = zeros(ncol,1);
  I = eye(ncol);
  converged = false;
  maxIter = inf;
  iter = 0;
  % alternatively update
  DDt = inv(D'*D+lambda1*I)*D';
  while ~converged
    iter = iter + 1;
    xtemp = x;
    x = DDt*(z-e);
    % x = (D'*D + lambda1*I)\(D'*(z-e));
    etemp = e;
    e = thres(z-D*x,lambda2);
    stopc = max(norm(e-etemp), norm(x-xtemp))/ndim;
    if stopc < 1e-6 || iter > maxIter
      converged = true;
    end
    % fval = func_proj(z,D,lambda1,lambda2,x,e);
    % fprintf('fval = %f\n', fval);
  end
end

function x = thres(y,mu)
  x = max(y-mu, 0);
  x = x + min(y + mu, 0);
end

function fval = func_proj(z,D,lambda1,lambda2,x,e)
  fval = 0;
  fval = fval + 0.5*norm(z-D*x-e)^2;
  fval = fval + 0.5*lambda1*norm(x)^2;
  fval = fval + lambda2*sum(abs(e));
end