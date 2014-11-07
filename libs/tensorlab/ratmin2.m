function [xy,v,g] = ratmin2(p,q,options)
%RATMIN2 Minimize bivariate and real polyanalytic rational functions.
%   [xy,v,g] = ratmin2(p,q) computes stationary points xy of a rational
%   function p/q in the two real variables (x,y) or in the complex variable
%   z (and conj(z)). The matrices p and q represent one of the following
%   two real polynomial types:
%
%      [1 y  y^2  ... y^dy]*p*[1;x;x^2;...;x^dx]      (analytic  bivariate)
%      [1 z' z'^2 ... z'^d]*p*[1;z;z^2;...;z^d]   (polyanalytic univariate)
%
%   where p(i,j) is the coefficient of the term x^j*y^i and z^j*conj(z)^i,
%   respectively. Each row of xy(i,:) is a stationary point of p/q, v(i) is
%   its function value p(xy(i,:))/q(xy(i,:)) and g(i,:) is the normalized
%   gradient
%
%      [g(xy(i,:))/|g|(|xy(i,:)|) h(xy(i,:))/|h|(|xy(i,:)|)],
%
%   where g := q^2*d(p/q)/dx or q^2*d(p/q)/dz and h := q^2*d(p/q)/dy or
%   q^2*d(p/q)/d(conj(z)) for bivariate and univariate p and q,
%   respectively.
%
%   ratmin2(p,q,options) may be used to set the following options:
%
%      options.<...>          - Parameters passed to polysol2, which solves
%                               the system q^2*d(p/q)/dx = q^2*d(p/q)/dy =
%                               0 or q^2*d(p/q)/dz = q^2*d(p/q)/d(conj(z))
%                               = 0 to obtain the solutions xy.
%      options.Univariate =   - True if p and q are to be interpreted as
%      p==p' && ~isreal(p) &&   polyanalytic univariate polynomials. False 
%      q==q' && ~isreal(q)      if they are bivariate polynomials.
%
%   See also ratmin, polyval2.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Numerical solution of
%       bivariate and polyanalytic polynomial systems", ESAT-SISTA Internal
%       Report 13-84, KU Leuven, 2013.

% Check the options structure.
pIsHermitian = size(p,1) == size(p,2) && all(all(tril(p)' == triu(p)));
qIsHermitian = size(p,1) == size(p,2) && all(all(tril(p)' == triu(p)));
if nargin < 3, options = struct; end
if ~isfield(options,'Univariate')
    options.Univariate = ~isreal(p) && ~isreal(q) && ...
                         pIsHermitian && qIsHermitian;
end

% Check the polynomials p and q.
if (~isreal(p) && ~pIsHermitian) || (~isreal(q) && ~qIsHermitian)
    error('ratmin2:pq','The polynomials p and q must be real-valued.');
end
if options.Univariate && ~(pIsHermitian && qIsHermitian)
    error('ratmin2:pq','The matrices p and q must be Hermitian.');
end

% If p and q are vectors, tell user to call ratmin.
if any(size(p) == numel(p)) && any(size(q) == numel(q))
    error('ratmin2:pq','Use ratmin for analytic univariate rationals.');
end

% Solve the system q^2*d(p/q)/dx = q^2*d(p/q)/dy = 0 or
% q^2*d(p/q)/dz = q^2*d(p/q)/dconj(z) = 0.
dpdx = p(:,2:end)*diag(1:size(p,2)-1);
dpdy = diag(1:size(p,1)-1)*p(2:end,:);
dqdx = q(:,2:end)*diag(1:size(q,2)-1);
dqdy = diag(1:size(q,1)-1)*q(2:end,:);
dpqdx = 0;
if ~isempty(dpdx), dpqdx = conv2(dpdx,q); end
if ~isempty(dqdx), dpqdx = dpqdx-conv2(p,dqdx); end
if size(p,2) == size(q,2), dpqdx = dpqdx(:,1:end-1); end
if options.Univariate
    dpqdy = dpqdx';
else
    dpqdy = 0;
    if ~isempty(dpdy), dpqdy = conv2(dpdy,q); end
    if ~isempty(dqdy), dpqdy = dpqdy-conv2(p,dqdy); end
    if size(p,1) == size(q,1), dpqdy = dpqdy(1:end-1,:); end
end
if nargout < 3
    xy = polysol2(dpqdx,dpqdy,options);
else
    [xy,g] = polysol2(dpqdx,dpqdy,options);
end

% Optionally compute the function values.
if nargout >= 2, v = polyval2(p,xy)./polyval2(q,xy); end

% Update plot legend.
if isfield(options,'Plot') && options.Plot
    if options.Univariate
        legend('q^2*d(p/q)/dz = 0','q^2*d(p/q)/dconj(z) = 0','z*');
    else
        legend('q^2*d(p/q)/dx = 0','q^2*d(p/q)/dy = 0','(x*,y*)');
    end
end
