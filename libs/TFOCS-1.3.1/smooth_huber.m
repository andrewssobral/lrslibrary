function op = smooth_huber(tau )

%SMOOTH_HUBER   Huber function generation.
%   FUNC = SMOOTH_QUAD( TAU ) returns a function handle that implements
%
%        FUNC(X) = 0.5 *( x.^2 )/tau               if |x| <= tau
%                = |x| - tau/2                     if |x| >  tau
%
%   All arguments are optional; the default value is tau = 1.
%   The Huber function has continuous gradient and is convex.
%
%   The function acts component-wise.  TAU may be either a scalar
%   or a vector/matrix of the same size as X
%
%   Does not support nonsmooth usage yet

% Does not yet fully support tfocs_tuples

if nargin == 0,
    tau = 1;
end
% op = @(varargin) smooth_huber_impl(tau, varargin{:} ); % old method

op = tfocs_smooth( @smooth_huber_impl );


% function [ v, g ] = smooth_huber_impl(tau, x, t ) % old method
function [ v, g ] = smooth_huber_impl(x)
  if nargin == 3,
      error( 'Proximity minimization not supported by this function.' );
  end
  if ~isscalar(tau) && ~size(tau) == size(x)
      error('smooth_huber: tau must be a scalar or the same size as the variable');
  end
  smallSet    = ( abs(x) <= tau );
  v           = smallSet.*( 0.5*(x.^2)./tau ) + (~smallSet).*( abs(x) - tau/2) ;
  if nargout > 1
      g   = sign(x).*min( 1, abs(x)./tau );
  end
end % new method

end % new method


% TFOCS v1.3 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2013 California Institute of Technology and CVX Research.
% See the file LICENSE for full license information.
