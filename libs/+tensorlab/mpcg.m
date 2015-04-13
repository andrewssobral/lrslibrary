function [x,flag,relres,iter] = mpcg(A,b,tol,maxit,M,~,x0)
%MPCG Modified preconditioned conjugate gradients method.
%   x = mpcg(A,b) attempts to solve the system of linear equations A*x = b
%   for x. The n-by-n coefficient matrix A must be symmetric and positive
%   definite. The column vector b must have length n. A can be a function
%   handle afun such that afun(x) returns A*x. MATLAB's pcg.m can in some
%   cases return x0 (which is often equal to 0) if it is the iterate with
%   the smallest residual. The only difference between mpcg and MATLAB's
%   pcg method is that this implementation returns the last iterate,
%   regardless of its residual.
%
%   mpcg(A,b,tol) specifies the tolerance of the method. If tol is [], then
%   mpcg uses the default, 1e-6.
%
%   mpcg(A,b,tol,maxit) specifies the maximum number of iterations. If
%   maxit is [], then mpcg uses the default, min(n,20).
%
%   mpcg(A,b,tol,maxit,M) uses a symmetric positive definite preconditioner
%   M and effectively solve the system inv(M)*A*x = inv(M)*b for x. If M is
%   [] then mpcg applies no preconditioner. M can be a function handle mfun
%   such that mfun(x) returns M\x. If M is 'SSOR', then mpcg applies a
%   Symmetric Successive Over-Relaxation preconditioner.
%
%   mpcg(A,b,tol,maxit,M,[],x0) specifies the initial guess. If x0 is [],
%   then mpcg uses the default, an all-zero vector.
%
%   [x,flag,relres,iter] = mpcg(A,b,...) also returns a convergence flag,
%   which is 0 if mpcg converged to the desired tolerance within maxit
%   iterations and 1 otherwise, the relative residual norm(A*x-b)/norm(b)
%   and the number of iterations, respectively.
%
%   See also bicgstab, cgs, gmres, lsqsr, qmr.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] J. Nocedal and S. J. Wright, "Numerical Optimization," Springer
%       Series in Operations Research, Springer, second ed., 2006.

% Check the options.
if nargin < 3 || isempty(tol), tol = 1e-6; end
if nargin < 4 || isempty(maxit), maxit = 20; end
maxit = min(maxit,length(b));
if nargin == 5 && ischar(M) && strcmpi(M,'SSOR') && isnumeric(A)
    M = @(x)tril(A)*((triu(A)*x)./diag(A));
end
PC = nargin > 4 && (isa(M,'function_handle') || ...
     (isnumeric(M) && all(size(M) == length(b))));

% Initialize PCG.
if nargin < 7 || isempty(x0)
    x = zeros(size(b));
    r = -b;
else
    x = x0;
    if isnumeric(A), r = A*x-b;
    else r = A(x)-b; end
end
if PC
    if isnumeric(M), y = M\r;
    else y = M(r); end
    d = -y;
    rr = r'*y;
else
    d = -r;
    rr = r'*r;
end
if nargin < 7 && ~PC, normb = sqrt(rr);
else normb = norm(b); end
flag = 1;

% Preconditioned conjugate gradients.
for iter = 1:maxit

    if isnumeric(A), Ad = A*d;
    else Ad = A(d); end
    alpha = rr/(d'*Ad);

    x = x+alpha*d;
    r = r+alpha*Ad;
    rr1 = rr;
    if PC
        if isnumeric(M), y = M\r;
        else y = M(r); end
        rr = r'*y;
    else
        rr = r'*r;
    end
    
    if PC, relres = norm(r)/normb;
    else relres = sqrt(rr)/normb; end
    if relres < tol, flag = 0; break; end

    beta = rr/rr1;
    if PC, d = -y+beta*d;
    else d = -r+beta*d; end

end
