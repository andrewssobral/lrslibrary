function [A_dual,E_dual,O_dual,iter,Y] = ThreeWayDec(D,MM, tau,lambda,tol, maxIter)
[m,n] = size(D);
DISPLAY_EVERY = 1;

% initialize
Y = D;
norm_two = norm(Y, 2);
Y = Y / norm_two;
A_dual = zeros(m, n);
E_dual = zeros(m, n);
O_dual = zeros(m, n);
mu = 1.25/norm(D) ;
rho = 1.25;
d_norm = norm(D, 'fro');
iter = 0;
converged = false;

% Invert confidence (object vs turbulence)
MMinv = 1 - MM;
while ~converged
  iter = iter + 1;
  
  temp_T = D - E_dual - O_dual + (1/mu)*Y;
  [U,S,V] = svd(temp_T, 'econ');
  diagS = diag(S);
  A_dual = U * diag(pos(diagS-1/mu)) * V';
  temp_T = D - A_dual - E_dual + (1/mu)*Y;
  O_dual = sign(temp_T) .* pos( abs(temp_T) - ((tau/mu).*MMinv) );
  temp_T = D - A_dual - O_dual + (1/mu)*Y;
  E_dual = ( 1/(1 + 2.*lambda/mu) ).*temp_T;
  Z = D - A_dual - E_dual - O_dual;
  Y = Y + mu*Z;
  obj_v = D(:)'*Y(:);
  mu = mu*rho;
  stoppingCriterion = norm(Z, 'fro') / d_norm;
  
  if mod( iter, DISPLAY_EVERY) == 0
    disp(['#Iteration ' num2str(iter) '  rank(A) ' num2str(rank(A_dual)) ...
      ' ||E||_0 ' num2str(length(find(abs(E_dual)>0)))...
      ' objvalue ' num2str(obj_v) '  Stopping Criterion ' ...
      num2str(stoppingCriterion)]);
  end
  if stoppingCriterion <= tol
    disp('3WD converged at:');
    disp(['Iteration ' num2str(iter) '  rank(A) ' num2str(rank(A_dual)) ...
      ' ||E||_0 ' num2str(length(find(abs(E_dual)>0)))  '  Stopping Criterion ' ...
      num2str(stoppingCriterion)]) ;
    converged = true ;
  end
  if ~converged && iter >= maxIter
    disp('Maximum iterations reached') ;
    converged = 1 ;
  end
end
