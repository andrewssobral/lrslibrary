function op = proj_l2( q )

%PROJ_L2   Projection onto the scaled 2-norm ball.
%    OP = PROJ_L2( Q ) returns an operator implementing the 
%    indicator function for the 2-norm ball of size q,
%    { X | norm( X, 2 ) <= q }. Q is optional; if omitted,
%    Q=1 is assumed. But if Q is supplied, it must be a positive
%    real scalar.
% Dual: prox_l2.m
% See also: prox_l2.m

if nargin == 0,
	q = 1;
elseif ~isnumeric( q ) || ~isreal( q )
	error( 'Argument must be a real scalar.' );
end
if numel(q) == 1
    if q <= 0, error('Argument must be positive'); end
    op = @(varargin)proj_l2_q( q, varargin{:} );
else
    if any( abs(q) < 10*eps ), error('Weight "q" must be nonzero'); end
    warning('TFOCS:experimental','Using experimental feature of TFOCS');
    
    op = @(varargin)proj_l2_qVec( q, varargin{:} );
end

function [ v, x ] = proj_l2_q( q, x, t )
v = 0;
switch nargin,
	case 2,
		if nargout == 2,
			error( 'This function is not differentiable.' );
		elseif norm( x(:), 'fro' ) > q, % GKC fix 2013 (for > 2D arrays)
			v = Inf;
		end
	case 3,
        nrm = norm(x(:),'fro'); % fixing, Feb '11, and GKC fix 2013
        if nrm > q
            x = x .* ( q / nrm );
        end
	otherwise,
		error( 'Not enough arguments.' );
end

% -- experimental version for when q is a vector --
function [ v, x ] = proj_l2_qVec( q, x, t )
v = 0;
switch nargin,
	case 2,
		if nargout == 2,
			error( 'This function is not differentiable.' );
		elseif norm( x./q, 'fro' ) > 1,
			v = Inf;
		end
	case 3,
        nrm = norm(x./q,'fro');
        if nrm > 1

            % We know x is of the form x0./( 1 + lambda*D2 )
            %   for some lambda > 0, but we don't have an easy
            %   way to know what lambda is.  So, treat this as
            %   a 1D minimization problem to find lambda.
            D = 1./(q);
            D2 = D.^2;
            Dx = D.*x;
            
%             lMax  = max( abs(x./D2) )*sqrt(numel(x));
            lMax  = 1.2*norm( abs(x./D2),'fro'); % a tighter bound
            fmin_opts  = optimset( 'TolX', 1e-10 );
%             MaxFunEvals: 500
%                    MaxIter:
            [lOpt,val,exitflag,output]    = fminbnd( @(l) (norm(Dx./(1+l*D2),'fro')-1)^2, 0, lMax,fmin_opts);
            if val > 1e-3, error('Proj_l2 failed to converge'); end
            x       = x./( 1 + lOpt*D2 );
        end
        
	otherwise,
		error( 'Not enough arguments.' );
end

% TFOCS v1.3 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2013 California Institute of Technology and CVX Research.
% See the file LICENSE for full license information.
