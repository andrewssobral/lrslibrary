function [U, Theta, V, numiter ] = OR1MP(m, n, r, Known, data, opts )
%ORTHOGONAL-RANK-1-MATRIX-PURSUIT Infinit dimension matching pursuit for 
%low rank matrix
%
%For the problem           min    L(X) = ||X - Y||^2
%                          s.t.   rank(X) <= r
%                          X = Theta * M = sum_i theta_i * M_i = U Theta V'
%   Detail variable
%   Input:
%         m ---- row number
%         n ---- column number
%         r ---- number of basis
%         Known ---- index of the known spot in the matrix
%         data ---- the content of the known spot in the matrix 
%         opts ---- parameters for the algorithm
%   Output:
%         U ---- output matrix: U
%         Theta ---- output matrix: Theta
%         V ---- output matrix: V
%         numiter ---- number of iterations
%
%   Copyright Zheng Wang @ Arizona State University
%   $Date: 2013/01/29$

fn = mfilename;
error(nargchk(5, 6, nargin));
if nargin < 6
    epsilon = 1e-4; % convergence threshold
else
    epsilon = opts.epsilon;
end

%addpath('PROPACK/');
%addpath('largescale_ops/');
%addpath('SLEP_package_4.1/');

% initialization, 
[indm, indn] = ind2sub([m, n], Known);
data( data == 0 )= eps;
res = sparse(indm, indn, data, m, n);
[indm, indn, data] = find(res);
U = [];
V = [];
Msup = [];

verbosity = 1;
printsyb = ['-', 'X', '|'];
if verbosity == 1
    fprintf('\nIteration:        ');
end

i = 0;
W = 0;
oldresnorm = 0;
gresnorm = 1;
yy = [];
nnorm = norm(data, 'fro');
% main iteration
% In OR1MP, the stop criterion is small gradient of residual
while (i < r) && (gresnorm > epsilon )
    % 1. find the top singular pair of the residual and update the gresnorm
    resvec = data - W;
    sparse_update(res, resvec); % sparse update the res using resvec

    [u, ~, v] = topsvd(res, 20); % run our power method for 10 iterations
    %[u, ~, v] = topsvd(res, 1);
    %[u, ~, v] = lansvd(res, 1, 'L'); % fast sparse svd using PROPACK
    %[u, s, v] = svds(res, 1); % use matlab sparse top svd
    
    resnorm = normest(res, 'fro')/nnorm;
    gresnorm = abs(resnorm - oldresnorm);
    oldresnorm = resnorm;

    % 2. update the weight Theta, the pursuit basis is uv', its weight is s.
    Mi = sparse_inp(u', v', indm, indn)';

    % b) use incremental inverse to solve the least sqare problem
    if i~=0
        Minv = inverse_incremental(Minv, Msup'*Mi, Mi'*Mi);
    else
        Minv = 1/(Mi'*Mi);
    end
    yy = [yy; Mi'*data];
    Theta = Minv*yy;
    Msup = [Msup Mi];
     
    U = [U u];
    V = [V v];

    % 3. update the learned matrix W = U' * diag(Theta) * V;
    W = Msup * Theta;
    
    if verbosity == 1
        fprintf('\b\b\b\b\b\b\b\b  %c:%4d', printsyb(1+mod(i,length(printsyb))), i);
    end
    i = i + 1;
end
% V = diag(Theta)*V;
numiter = i;
% fprintf( '\n OR1MP run %d rounds! \n', numiter);

% %* Sparse selection: for sharpe low rank problem, we may need the lasso 
% %fine selection. Use lasso to learn a more sparse weights.
% % Starting point
% opts.init=2;        % starting from a zero point
% % termination criterion
% opts.tFlag=5;       % run .maxIter iterations
% opts.maxIter=100;   % maximum number of iterations
% % normalization
% opts.nFlag=0;       % without normalization
% % regularization
% opts.rFlag=1;       % the input parameter 'rho' is a ratio in (0, 1)
% opts.mFlag=0;       % treating it as compositive function
% opts.lFlag=0;       % Nemirovski's line search
% % get the final sparse Theta by lasso in SLEP toolbox
% [Theta, ~, ~] = LeastR(Msup, data, 0.00001, opts); %  [W, funVal, ValueL] 00001

function Ninv = inverse_incremental(Minv, MMi, d)
% calculate the inverse of the blocked matrix of 
P = MMi' * Minv; % vector
q = 1/(d - P*MMi); % scaler
y = q*P;
Ninv = [Minv+P'*y, -y'; -y, q];

function [u, s, v] = topsvd(A, round)
% calculate the top svd of matrix A using round iterations
stopeps = 1e-3;
[m,n] = size(A);
u = ones(m,1); % this is the sigma*u
vo = 0;
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
u = u/nu;
v = v'/nv;
% v = v/nv;
s = nu*nv;
