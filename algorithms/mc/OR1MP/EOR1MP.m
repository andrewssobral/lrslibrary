function [U, Theta, V, numiter ] = EOR1MP(m, n, r, Known, data, opts )
%MP-MATRIX-COMPLETION Infinite dimension MP for matrix completion 
%
%For the problem           min  sum L(X),
%                          s.t. X = Theta * M = sum_i theta_i * M_i.
%   Detailed explanation goes here
%   Input:
%         m ---- row number
%         n ---- column number
%         r ---- number of basis
%         Known ---- index of the known elements in the matrix
%         data ---- the content of the known elements in the matrix 
%         opts ---- parameters for the algorithm
%   Output:
%         U ---- output matrix: U
%         V ---- output matrix: V'
%         Theta ---- output weights 
%         numiter ---- number of iterations
%
%   Copyright: Zheng Wang @ ASU
%   $Date: 2013/12/14$

% addpath('largescale_ops/');
fn = mfilename;
error(nargchk(5, 6, nargin));
if nargin < 6
    verbosity = 1;
else
    verbosity = opts.verbosity;
end
[indm, indn] = ind2sub([m, n], Known);
data(data == 0)= eps;
res = sparse(indm, indn, data, m, n);
[indm, indn, data] = find(res);

% main iteration
U     = [];
V     = [];
Msup  = [];
i     = 0;
W     = 0;
printsyb  = ['-', 'X', '|'];
if verbosity == 1
    fprintf('\nIteration:        ');
end

% In regerssion pursuit, the stop criterion is small gradient of residual
Theta = [];
while (i < r)
    % 1. find the top singular pair of the residual
    resvec = data - W; sparse_update(res, resvec);    
    [u, ~, v] = topsvd(res, 20);

    % 2. update the weight Theta, the pursuit basis is uv', its weight is s.
    Mi = sparse_inp(u', v', indm, indn)';
    U    = [U u];
    V    = [V v];
    Msup = [Msup Mi];
    Sol  = inv(Msup'*Msup) * Msup' * data; % optimal line search by least square   
    if i == 0
        Theta = Sol;
    else
        Theta = [Theta * Sol(1) Sol(2)];
    end    
    W    = Msup * Sol; 
    Msup = W;
    if verbosity == 1
        fprintf('\b\b\b\b\b\b\b\b  %c:%4d', printsyb(1+mod(i,length(printsyb))), i);
    end
    i = i + 1;
end
numiter = i;
fprintf( '\n EOR1MP run %d rounds! \n', numiter);

function [u, s, v] = topsvd(A, round)
% calculate the top svd of matrix A using round iterations
stopeps = 1e-3;
[m, n]  = size(A);
u       = ones(m,1);
vo      = 0;
for i=1:round
    v = u'*A/(norm(u))^2;
    u = A*v'/(norm(v))^2;
    if norm(v-vo) < stopeps
        break
    end
    vo = v;
end
nu = norm(u);
nv = norm(v);
u  = u/nu;
v  = v'/nv;
s  = nu*nv;
