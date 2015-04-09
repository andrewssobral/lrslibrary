function Hv = hessvec_fd(v,FUN,x,gx,s)
%HESSVEC_FD   Hessian vector product finite difference approximation.
%
%   HV = HESSVEC_FD(V,FUN,X) computes a forward finite difference
%   approximation of the Hessian vector product, H(X)*V, where V is given 
%   and H(X) is the Hessian of the function FUN at the point X. The
%   approximation is given by
%
%              G(X+S*V) - G(X)
%     H(X)*V = ---------------
%                     S
%
%   where G(X) is the gradient of the function FUN at the point X and 
%   h is the difference step (default = 1e-8*(1+norm(X)).
%
%   HV = HESSVEC(V,FUN,X,GX) uses GX for the value of G(X), in the
%   case that is has already been computed outside of this method.
%
%   HV = HESSVEC(V,FUN,X,GX,S) uses S as the difference step.
%
%   This method should not be called directly, but only by Poblano Toolbox
%   optimization algorithms.
%
%   The number of calls to FUN is kept track of in the global variable
%   NFEV_HESSVEC_FD to facilitate using this method as a callback function
%   (e.g., from an iterative linear solver such as Matlab's SYMMLQ).
%
%   See also TN.
%
%MATLAB Poblano Toolbox.
%Copyright 2009-2012, Sandia Corporation.

%% Number of calls to FUN
global nfev_hessvec_fd;

%% Check input
if nargin < 3
    error('HESSVEC_FD => at least two input arguments are required (V,FUN,X)');
end

% Compute gradient if not given
if nargin < 4
    [f,gx] = feval(FUN,x);
    nfev_hessvec_fd = nfev_hessvec_fd + 1;
end

% Compute a difference step if not given
if nargin < 5
    s = 1e-8*(1+norm(x));
end

% Compute the gradient at the new point
[f,gxsv] = feval(FUN,x+s*v);
nfev_hessvec_fd = nfev_hessvec_fd + 1;

%% Hessian vector product approximation
Hv = (gxsv-gx)/s;
