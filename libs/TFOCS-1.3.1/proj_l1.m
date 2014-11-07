function op = proj_l1( q )
%PROJ_L1   Projection onto the scaled 1-norm ball.
%    OP = PROJ_L1( Q ) returns an operator implementing the 
%    indicator function for the 1-norm ball of radius q,
%    { X | norm( X, 1 ) <= q }. Q is optional; if omitted,
%    Q=1 is assumed. But if Q is supplied, it must be a positive
%    real scalar.
% Dual: prox_linf.m
% See also: prox_linf, prox_l1, proj_linf

if nargin == 0,
	q = 1;
elseif ~isnumeric( q ) || ~isreal( q ) || numel( q ) ~= 1 || q <= 0,
	error( 'Argument must be positive.' );
end
op = @(varargin)proj_l1_q(q, varargin{:} );

function [ v, x ] = proj_l1_q( q, x, t )
v = 0;
switch nargin,
case 2,
	if nargout == 2,
		error( 'This function is not differentiable.'  );
	elseif norm( x(:), 1 ) > q,
		v = Inf;
	end
case 3,
    s      = sort(abs(nonzeros(x)),'descend');
    cs     = cumsum(s);
    % ndx    = find( cs - (1:numel(s))' .* [ s(2:end) ; 0 ] >= q, 1 );
    ndx    = find( cs - (1:numel(s))' .* [ s(2:end) ; 0 ] >= q+2*eps(q), 1 ); % For stability
    if ~isempty( ndx )
        thresh = ( cs(ndx) - q ) / ndx;
        x      = x .* ( 1 - thresh ./ max( abs(x), thresh ) ); % May divide very small numbers
    end
otherwise,
    error( 'Not enough arguments.' );
end

% TFOCS v1.3 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2013 California Institute of Technology and CVX Research.
% See the file LICENSE for full license information.
