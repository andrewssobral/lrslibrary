function [x,v,g] = ratmin(p,q,options)
%RATMIN Minimize a rational function.
%   [x,v,g] = ratmin(p,q) computes stationary points x of a real rational
%   function p/q. The vectors p and q represent the polynomials
%   p*[x^dp;...;1] and q*[x^dq;...;1], respectively. I.e., p(d-i+1) and
%   q(d-i+1) are the coefficients of the term x^i in p and q, respectively.
%   Each x(i) is a stationary point of p/q, v(i) is its function value
%   p(x(i))/q(x(i)) and g(i) is the normalized gradient g(x(i))/|g|(|x(i)|)
%   where g is q^2*d(p/q)/dx.
%
%   ratmin(p,q,options) may be used to set the following options:
%
%      options.TolReal = 1e2 - A root x* of q^2*d(p/q)/dx is accepted as a
%                              stationary point if its imaginary part
%                              is smaller than options.TolReal* ...
%                              eps(real(x*)).
%
%   See also ratmin2, polymin, polyval.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Check the options structure.
if nargin < 3, options = struct; end
if ~isfield(options,'TolReal'), options.TolReal = 1e2; end

% Compute and normalize the gradient q^2*d(p/q)/dx.
p = p(:).'; q = q(:).';
dpdx = polyder(p);
dqdx = polyder(q);
g = conv(dpdx,q)-conv(p,dqdx);
if length(p) == length(q), g = g(2:end); end
if ~any(isinf(g/max(abs(g)))), g = g/max(abs(g)); end

% Compute the stationary points and filter complex solutions.
x = roots(g);
idx = abs(imag(x)) < options.TolReal*eps(real(x));
x = real(x(idx));

% Optionally compute the function value and gradient.
if nargout >= 2, v = polyval(p,x)./polyval(q,x); end
if nargout == 3, g = polyval(g,x)./polyval(abs(g),abs(x)); end
