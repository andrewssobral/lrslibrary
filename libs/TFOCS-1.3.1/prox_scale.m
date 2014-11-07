function op = prox_scale( proxF, scale )

%PROX_SCALE   Scaling a proximity/projection function.
%    PSCALE = PROX_SCALE( PROXF, s ) is the proximity function formed
%    by multiplying the input by the real value SCALE, then calling the
%    function PROXF. In other words,
%        PSCALE( y ) = PROXF( s * y ).
%    In three-argument form, [v,x]=PSCALE(y,t) finds the minimizer of
%        PROXF( s * x ) + (0.5/t)*||x-y||_2^2
%
%    Why scale the input, and not the output? When the proximity function
%    is a norm, the two choices are equivalent when SCALE >= 0, since
%    ||s*x|| = abs(s)*||x||. However, for *indicator* functions, scaling
%    the output is useless. Scaling the input grows or shrinks the domain
%    of PSCALE by choosing SCALE<1 or SCALE>1, respectively.
%
%    Note that when constructing our equivalence we assume a particular
%    form for the proximity minimization. More generally, PROJ_SCALE will
%    preserve equivalence if 0.5*||x-y||^2 is replaced with D(x-y), as
%    long as D() satisfies the following homogeneity property:
%        D(s*z) = s^2*D(z)
%    If this is not the case, you will not be able to use the TFOCS
%    internal scaling with this function.

if ~isa( proxF, 'function_handle' ),
    error( 'The first argument must be a function handle.' );
elseif ~isa( scale, 'double' ) || ~isreal( scale ) || numel( scale ) ~= 1 || scale == 0,
    error( 'The second argument must be a nonzero real scalar.' );
elseif scale == 1,
	op = proxF;
else
	op = @(varargin)prox_scale_impl( proxF, scale, varargin{:} );
end

function varargout = prox_scale_impl( proxF, s, y, t )
no = max(nargout,1);
if nargin < 4 || t == 0,
    [ varargout{1:no} ] = proxF( s * y );
else
    % For three-input mode, we wish to simulate
    %   x = argmin_x projectorF(s*x) + (0.5/t)*||x-y||^2
    % Setting z = s * x, we have
    %   z = argmin_z projectorF(z) + (0.5/(s^2*t))*||z-s*y||^2
    % Therefore, we have
    %   [v,z] = projectorF( s * y, t * s^2 ); x = z / s;
	[ varargout{1:no} ] = proxF( s * y, s^2 * t );
end
if no > 1,
    varargout{2} = varargout{2} / s;
end

% TFOCS v1.3 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2013 California Institute of Technology and CVX Research.
% See the file LICENSE for full license information.
