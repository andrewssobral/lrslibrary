function [xy,v,g] = polymin2(p,options)
%POLYMIN2 Minimize bivariate and real polyanalytic polynomials.
%   [xy,v,g] = polymin2(p) computes stationary points xy of a real
%   polynomial p in the two real variables (x,y) or in the complex variable
%   z (and conj(z)). The matrix p represents one of the following two real
%   polynomial types:
%
%      [1 y  y^2  ... y^dy]*p*[1;x;x^2;...;x^dx]      (analytic  bivariate)
%      [1 z' z'^2 ... z'^d]*p*[1;z;z^2;...;z^d]   (polyanalytic univariate)
%
%   where p(i,j) is the coefficient of the term x^j*y^i and z^j*conj(z)^i,
%   respectively. Each row of xy(i,:) is a stationary point of p, v(i) is
%   its function value p(xy(i,:)) and g(i,:) is the normalized gradient
%
%      [g(xy(i,:))/|g|(|xy(i,:)|) h(xy(i,:))/|h|(|xy(i,:)|)],
%
%   where g := dp/dx or dp/dz and h := dp/dy or dp/d(conj(z)) for bivariate
%   and univariate p, respectively.
%
%   polymin2(p,options) may be used to set the following options:
%
%      options.<...>        - Parameters passed to polysol2, which solves
%                             the system dp/dx = dp/dy = 0 or dp/dz =
%                             dp/d(conj(z)) = 0 to obtain the solutions xy.
%      options.Univariate = - True if p is to be interpreted as a
%      p==p' && ~isreal(p)    polyanalytic univariate polynomial. False if
%                             it is a bivariate polynomial.
%
%   See also polymin, polyval2.

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
if nargin < 2, options = struct; end
if ~isfield(options,'Univariate')
    options.Univariate = ~isreal(p) && pIsHermitian;
end

% Check the polynomial p.
if ~isreal(p) && ~pIsHermitian
    error('polymin2:pq','The polynomial p must be real-valued.');
end
if options.Univariate && ~pIsHermitian
    error('polymin2:pq','The matrix p must be Hermitian.');
end

% If p is a vector, tell user to call polymin.
if any(size(p) == numel(p))
    error('polymin2:p','Use polymin for analytic univariate polynomials.');
end

% Solve the system dp/dx = dp/dy = 0 or dp/dz = dp/dconj(z) = 0.
dpdx = p(:,2:end)*diag(1:size(p,2)-1);
dpdy = diag(1:size(p,1)-1)*p(2:end,:);
if nargout < 3
    xy = polysol2(dpdx,dpdy,options);
else
    [xy,g] = polysol2(dpdx,dpdy,options);
end

% Optionally compute the function values.
if nargout >= 2, v = polyval2(p,xy); end

% Update plot legend.
if isfield(options,'Plot') && options.Plot
    if options.Univariate
        legend('dp/dz = 0','dp/dconj(z) = 0','z*');
    else
        legend('dp/dx = 0','dp/dy = 0','(x*,y*)');
    end
end
