function [x,v,g] = polymin(p,options)
%POLYMIN Minimize a polynomial.
%   [x,v,g] = polymin(p) computes stationary points x of a real polynomial
%   p, which is a vector that represents the polynomial p*[x^d;...;1],
%   i.e., p(d-i+1) is the coefficient of the term x^i. Each x(i) is a
%   stationary point of p, v(i) is its function value p(x(i)) and g(i) is
%   the normalized gradient g(x(i))/|g|(|x(i)|) where g is dp/dx.
%
%   polymin(p,options) may be used to set the following options:
%
%      options.TolReal = 1e2 - A root x* of dp/dx is accepted as a
%                              stationary point if its imaginary part
%                              is smaller than options.TolReal* ...
%                              eps(real(x*)).
%
%   See also polymin2, polyval.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Check the options structure.
if nargin < 2, options = struct; end
if ~isfield(options,'TolReal'), options.TolReal = 1e2; end

% Compute and normalize the gradient dp/dx.
p = p(:).';
g = polyder(p);
if ~any(isinf(g/max(abs(g)))), g = g/max(abs(g)); end

% Compute the stationary points and filter complex solutions.
x = roots(g);
idx = abs(imag(x)) < options.TolReal*eps(real(x));
x = real(x(idx));

% Optionally compute the function value and gradient.
if nargout >= 2, v = polyval(p,x); end
if nargout == 3, g = polyval(g,x)./polyval(abs(g),abs(x)); end
