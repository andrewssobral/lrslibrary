function v = polyval2(p,x,y)
%POLYVAL2 Evaluate bivariate and univariate polyanalytic polynomials.
%   v = polyval2(p,xy) computes the values v(i) of the polynomial p in two
%   real variables (x,y) at the points (xy(i,1),xy(i,2)). The polynomial is
%   represented as [1 y ... y^dy]*p*[1;x;...;x^dx].
%
%   v = polyval2(p,x,y) computes the values v(i,j) of the polynomial p in
%   two real variables (x,y) at the points (x(j),y(i)). The polynomial is
%   represented as [1 y ... y^dy]*p*[1;x;...;x^dx].
%
%   v = polyval2(p,z) computes the values v(i) of the polynomial p in the 
%   complex variable (z,conj(z)) at the points (z(i),conj(z(i))). The
%   polynomial is represented as [1 z' ... z'^d1]*p*[1;z;...;z^d2].
%
%   See also polyval, polymin2.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

if isempty(p), p = 0; end
if nargin == 2
    
    % Points supplied as [x y] or as x and implicitly y := conj(x).
    % Evaluate polynomial in (x(i),y(i)).
    if size(x,2) == 2
        y = x(:,2);
        x = x(:,1).';
        isReal = false;
    else
        x = x(:).';
        y = x';
        isReal = size(p,1) == size(p,2) && all(all(tril(p)' == triu(p))); 
    end
    q = p(:,end);
    if size(p,2) == 1, q = repmat(q,1,length(x)); end
    for j = size(p,2)-1:-1:1
        q = bsxfun(@plus,p(:,j),bsxfun(@times,q,x));
    end
    q = q.';
    v = q(:,end);
    for i = size(q,2)-1:-1:1
        v = q(:,i)+v.*y;
    end
    if isReal
        v = real(v);
    end
    
else
    
    % Points supplied as x and y. Evaluate polynomial in all (x(j),y(i)).
    x = x(:).';
    y = y(:).';
    q = p(:,end);
    if size(p,2) == 1, q = repmat(q,1,length(x)); end
    for j = size(p,2)-1:-1:1
        q = bsxfun(@plus,p(:,j),bsxfun(@times,q,x));
    end
    q = q.';
    v = q(:,end);
    if size(q,2) == 1, v = repmat(v,1,length(y)); end
    for i = size(q,2)-1:-1:1
        v = bsxfun(@plus,q(:,i),bsxfun(@times,v,y));
    end
    v = v.';
    
end
