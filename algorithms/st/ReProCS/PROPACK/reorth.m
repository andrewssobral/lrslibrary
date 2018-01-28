function [r,normr,nre,s] = reorth(Q,r,normr,index,alpha,method)

%REORTH   Reorthogonalize a vector using iterated Gram-Schmidt
%
%   [R_NEW,NORMR_NEW,NRE] = reorth(Q,R,NORMR,INDEX,ALPHA,METHOD)
%   reorthogonalizes R against the subset of columns of Q given by INDEX. 
%   If INDEX==[] then R is reorthogonalized all columns of Q.
%   If the result R_NEW has a small norm, i.e. if norm(R_NEW) < ALPHA*NORMR,
%   then a second reorthogonalization is performed. If the norm of R_NEW
%   is once more decreased by  more than a factor of ALPHA then R is 
%   numerically in span(Q(:,INDEX)) and a zero-vector is returned for R_NEW.
%
%   If method==0 then iterated modified Gram-Schmidt is used.
%   If method==1 then iterated classical Gram-Schmidt is used.
%
%   The default value for ALPHA is 0.5. 
%   NRE is the number of reorthogonalizations performed (1 or 2).

% References: 
%  Aake Bjorck, "Numerical Methods for Least Squares Problems",
%  SIAM, Philadelphia, 1996, pp. 68-69.
%
%  J.~W. Daniel, W.~B. Gragg, L. Kaufman and G.~W. Stewart, 
%  ``Reorthogonalization and Stable Algorithms Updating the
%  Gram-Schmidt QR Factorization'', Math. Comp.,  30 (1976), no.
%  136, pp. 772-795.
%
%  B. N. Parlett, ``The Symmetric Eigenvalue Problem'', 
%  Prentice-Hall, Englewood Cliffs, NJ, 1980. pp. 105-109

%  Rasmus Munk Larsen, DAIMI, 1998.

% Check input arguments.
% warning('PROPACK:NotUsingMex','Using slow matlab code for reorth.')
if nargin<2
  error('Not enough input arguments.')
end
[n k1] = size(Q);
if nargin<3 | isempty(normr)
%  normr = norm(r);
  normr = sqrt(r'*r);
end
if nargin<4 | isempty(index)
  k=k1;
  index = [1:k]';
  simple = 1;
else
  k = length(index);
  if k==k1 & index(:)==[1:k]'
    simple = 1;
  else
    simple = 0;
  end
end
if nargin<5 | isempty(alpha)
  alpha=0.5; % This choice garanties that 
             % || Q^T*r_new - e_{k+1} ||_2 <= 2*eps*||r_new||_2,
             % cf. Kahans ``twice is enough'' statement proved in 
             % Parletts book.
end
if nargin<6 | isempty(method)
   method = 0;
end
if k==0 | n==0
  return
end
if nargout>3
  s = zeros(k,1);
end


normr_old = 0;
nre = 0;
while normr < alpha*normr_old | nre==0
  if method==1
    if simple
      t = Q'*r;
      r = r - Q*t;
    else
      t = Q(:,index)'*r;
      r = r - Q(:,index)*t;
    end
  else    
    for i=index, 
      t = Q(:,i)'*r; 
      r = r - Q(:,i)*t;
    end
  end
  if nargout>3
    s = s + t;
  end
  normr_old = normr;
%  normr = norm(r);
  normr = sqrt(r'*r);
  nre = nre + 1;
  if nre > 4
    % r is in span(Q) to full accuracy => accept r = 0 as the new vector.
    r = zeros(n,1);
    normr = 0;
    return
  end
end
