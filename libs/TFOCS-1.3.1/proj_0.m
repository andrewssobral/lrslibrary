function op = proj_0(offset)

%PROJ_0     Projection onto the set {0}
%    OP = PROJ_0 returns an implementation of the indicator 
%    function for the set including only zero.
%
%    OP = PROJ_0( c ) returns an implementation of the
%    indicator function of the set {c}
%    If c is a scalar, this is interpreted as c*1
%    where "1" is the all ones object of the appropriate size.
%
% See also prox_0.m, proj_Rn.m and smooth_constant.m,
%   which are the Fenchel conjugates of this function.

if nargin == 0
    op = @proj_0_impl;
else
    op = @(varargin) proj_0_impl_q( offset, varargin{:} );
end

function [ v, x ] = proj_0_impl( x, t )
v = 0;
switch nargin,
	case 1,
		if nargout == 2,
			error( 'This function is not differentiable.' );
        elseif any( x(:) ),
            v = Inf;
        end
	case 2,
        % "t" variable has no effect
		x = 0*x;
	otherwise,
		error( 'Not enough arguments.' );
end
function [ v, x ] = proj_0_impl_q( c,  x, t )
v = 0;
switch nargin,
	case 2,
		if nargout == 2,
			error( 'This function is not differentiable.' );
        end
        if isscalar(c) 
            if any( x(:) - c )
                v = Inf;
            end
        elseif any( x(:)  - c(:) ),
            v = Inf;
        end
	case 3,
        % "t" variable has no effect
        if isscalar(c) && ~isscalar(x)
            x = c*ones( size(x) );
        else
            x = c;
        end

	otherwise,
		error( 'Not enough arguments.' );
end
% TFOCS v1.3 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2013 California Institute of Technology and CVX Research.
% See the file LICENSE for full license information.
