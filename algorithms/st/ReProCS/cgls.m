function [x,flag,resNE,iter] = cgls(A,b,shift,tol,maxit,prnt,x0)
%CGLS Conjugate Gradient Least Squares
%   X = CGLS(A,B) attempts to solve the system of linear equations A*X=B
%   for X. The M-by-N coefficient matrix A and right hand side column
%   vector B of length N are required input arguments.
%
%   X = CGLS(AFUN,B) accepts a function handle AFUN instead of the matrix
%   A. AFUN(X,1) accepts a vector input X and returns the matrix-vector
%   product A*X. AFUN(X,2) returns the matrix-vector product A'*X instead.
%   In all of the following syntaxes, you can replace A by AFUN.
%
%   X = CGLS(A,B,SHIFT) specifies a regularization parameter SHIFT. If
%   SHIFT is 0, then CGLS is Hestenes and Stiefel's specialized form of the
%   conjugate-gradient method for least-squares problems. If SHIFT is
%   nonzero, the system (A'*A + SHIFT*I)*X = A'*B is solved. Here I is the
%   N-by-N identity matrix.
%
%   X = CGLS(A,B,SHIFT,TOL) specifies the tolerance of the method. If TOL
%   is [] then CGLS uses the default, 1e-6.
%
%   X = CGLS(A,B,SHIFT,TOL,MAXIT) specifies the maximum number of
%   iterations. If MAXIT is [] then CGLS uses the default, 20.
%
%   X = CGLS(A,B,SHIFT,TOL,MAXIT,PRNT) specifies if output should be
%   generated during each iteration (PRNT == true). If PRNT is [] then no 
%   output is given.
%
%   X = CGLS(A,B,SHIFT,TOL,MAXIT,PRNT,X0) specifies the N-by-1 initial
%   solution that is used. If X0 is [] then CGLS uses the default, 
%   X0 = zeros(N,1).
%
%   [X,FLAG] = CGLS(A,B,...) also returns a convergence FLAG:
%    1. CGLS converged to the desired tolerance TOL within MAXIT
%       iterations.
%    2. CGLS iterated MAXIT times but did not converge.
%    3. Matrix (A'*A + SHIFT*I) seems to be singular or indefinite.
%    4. Instability seems likely meaning (A'*A + SHIFT*I) indefinite and 
%       NORM(X) decreased.
%
%   [X,FLAG,RESNE] = CGLS(A,B,...) also returns the relative residual for 
%   the normal equations NORM(A'*B - (A'*A + SHIFT*I)*X)/NORM(A'*B).
%
%   [X,FLAG,RESNE,ITER] = CGLS(A,B,...) also returns the iteration number 
%   at which X was computed: 0 <= ITER <= MAXIT.
%
%   See also LSQR, PCG, FUNCTION_HANDLE.

%   01 Sep 1999: First version.
%                Per Christian Hansen (DTU) and Michael Saunders (visiting
%                DTU).
%   22 Jan 2013: Updated syntax and documentation.
%                Folkert Bleichrodt (CWI).


% Assign default values to unspecified parameters
if (nargin < 3 || isempty(shift)), shift = 0;    end
if (nargin < 4 || isempty(tol))  , tol   = 1e-6; end
if (nargin < 5)                  , maxit = [];   end
if (nargin < 6 || isempty(prnt)) , prnt  = 0;    end
if (nargin < 7)                  , x0 = [];      end


if isa(A, 'numeric')
    explicitA = true;
elseif isa(A, 'function_handle')
    explicitA = false;
else
    error('A must be numeric or a function handle.');
end

% handle initial guess, if passed as argument
if explicitA
    [m,n] = size(A);

    if ~isempty(x0)
        x = x0;
    else
        x = zeros(n,1);
    end
    
    r = b - A*x;
    s = A'*r-shift*x;
else
    m = size(b,1);
    
    if ~isempty(x0)
        x = x0;
        r = b - A(x,1);
        s = A(r,2) - shift*x;
        n = size(s,1);
    else
        r = b;
        s = A(b,2);
        n = size(s,1);
        x = zeros(n,1);
    end
    
end

% determine default for maxit
if isempty(maxit)
    maxit = min([m,n,20]);
end

% Initialize
p      = s;
norms0 = norm(s);
gamma  = norms0^2;
normx  = norm(x);
xmax   = normx;
k      = 0;
flag   = 0;

if prnt
    head = '    k       x(1)             x(n)           normx        resNE';
    form = '%5.0f %16.10g %16.10g %9.2g %12.5g\n';
    disp('  ');   disp(head);
    fprintf(form, k, x(1), x(n), normx, 1);
end

indefinite = 0;

%--------------------------------------------------------------------------
% Main loop
%--------------------------------------------------------------------------
while (k < maxit) && (flag == 0)
    
    k = k+1;
    
    % q = A p
    if explicitA
        q = A*p;
    else
        q = A(p,1);
    end
        
    delta = norm(q)^2  +  shift*norm(p)^2;
    if delta <= 0, indefinite = 1;   end
    if delta == 0, delta      = eps; end
    alpha = gamma / delta;
    
    x     = x + alpha*p;
    r     = r - alpha*q;
       
    if explicitA
        s = A'*r - shift*x;
    else
        s = A(r,2) - shift*x;
    end
   
    norms  = norm(s);
    gamma1 = gamma;
    gamma  = norms^2;
    beta   = gamma / gamma1;
    p      = s + beta*p;
    
    % Convergence
    normx = norm(x);
    xmax  = max(xmax, normx);
    flag  = (norms <= norms0 * tol) || (normx * tol >= 1);
    
    % Output
    resNE = norms / norms0;
    if prnt, fprintf(form, k, x(1), x(n), normx, resNE); end
end % while

iter = k;

shrink = normx/xmax;
if k == maxit,          flag = 2; end
if indefinite,          flag = 3; end
if shrink <= sqrt(tol), flag = 4; end
