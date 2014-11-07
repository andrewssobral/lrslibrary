function [ x, odata, opts ] = solver_L1RLS( A, b, lambda, x0, opts )
% SOLVER_L1RLS l1-regularized least squares problem, sometimes called the LASSO.
% [ x, odata, opts ] = solver_L1RLS( A, b, lambda, x0, opts )
%    Solves the l1-regularized least squares problem
%        minimize (1/2)*norm( A * x - b )^2 + lambda * norm( x, 1 )
%    using the Auslender/Teboulle variant with restart. A must be a matrix
%    or a linear operator, b must be a vector, and lambda must be a real
%    positive scalar. The initial point x0 and option structure opts are
%    both optional.
%
%   Note: this formulation is sometimes referred to as "The Lasso"

error(nargchk(3,5,nargin));
if nargin < 4, x0 = []; end
if nargin < 5, opts = []; end
if ~isfield( opts, 'restart' ), 
    opts.restart = 100; 
end

[x,odata,opts] = tfocs( smooth_quad, { A, -b }, prox_l1( lambda ), x0, opts );

% TFOCS v1.3 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2013 California Institute of Technology and CVX Research.
% See the file LICENSE for full license information.

