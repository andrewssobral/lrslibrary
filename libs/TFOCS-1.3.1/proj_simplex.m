function op = proj_simplex( q )

%PROJ_SIMPLEX	Projection onto the simplex.
%    OP = PROJ_SIMPLEX( Q ) returns an nonsmooth function that
%    represents the scaled simplex { x | x >= 0, sum(x) <= q }.
%    Q is optional; if not supplied, it defaults to 1. If it  is
%    supplied, it must be a real positive scalar.
%
%   See also proj_psdUTrace.m (the matrix-analog of this function)

if nargin == 0,
	q = 1;
elseif ~isnumeric( q ) || ~isreal( q ) || numel( q ) ~= 1 || q <= 0,
	error( 'Argument must be positive.' );
end
op = @(varargin)proj_simplex_q( q,varargin{:} );

function [ v, x ] = proj_simplex_q( q, x, t )
v = 0;

% We have two options when input x is a matrix:
%   Does the user want to treat it as x(:), i.e. "vectorize" it ?
%   Or does the user want to treat each column separately?
% Most other functions (e.g. l1, linf) treat it as x(:)
% so that will be the default. However, we leave it
% as a hard-coded option so that the user can change it
% if they want.
VECTORIZE = true;

if nargin > 2 && t > 0,
	if any( x(:) < 0 ) || any( sum( x ) > q ),
        if VECTORIZE
            s     = sort( x(:), 'descend' );
        else
            s     = sort( x, 'descend' );
        end
        if q < eps(s(1))
            error('Input is scaled so large compared to q that accurate computations are difficult');
            % since then cs(1) = s(1) - q is  not even guaranteed to be
            % smaller than q !
        end
        if size(x,2) == 1 || VECTORIZE
            cs    = ( cumsum(s) - q ) ./ ( 1 : numel(s) )';
            ndx   = nnz( s > cs );
            x     = max( x - cs(ndx), 0 );
        else
            % vectorized for several columns
            cs    = diag( 1./( 1 : size(s,1) ) )*( cumsum(s) - q );
            ndx   = sum( s > cs );
            ndx   = sub2ind( size(x), ndx, 1:size(x,2)  );
            x     = max( x - repmat(cs(ndx),size(x,1),1), 0 );
        end
	end
% Sept 6 2012, adding factor of 100 in the line below:
elseif any( x(:) < 0 ) || any( abs( sum(x) / q - 1 ) > sqrt(numel(x)) * 100 *eps )
    v = Inf;
end

% TFOCS v1.3 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2013 California Institute of Technology and CVX Research.
% See the file LICENSE for full license information.
