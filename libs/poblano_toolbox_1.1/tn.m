function out = tn(FUN,x0,varargin)
%TN   Truncated Newton minimization.
%
%   OUT = TN(FUN,X0) minimizes FUN starting at the point X0 using a 
%   Hessian-free truncated Newton method. The outer loop is Newton's
%   method, and the inner loop uses a conjugate gradient method to compute
%   an approximation to the Newton direction. Furthermore, the Hessian
%   vector product is approximated using forward finite differences using 
%   HESSVEC_FD. FUN is a handle for a function that takes a single vector
%   input and returns two arguments --- the scalar function value and the 
%   vector-valued gradient. See POBLANO_OUT for details of the output
%   parameters.
%
%   OUT = TN(FUN,X0,'param',value,...) specifies a
%   parameters and its value. See POBLANO_PARAMS for further details on
%   standard parameters. Additionally, TN requires
%
%   'CGSolver' - Matlab CG method to use {'symmlq'}
%     'symmlq' : symmlq (designed for symmetric indefinite systems) 
%     'pcg' :    pcg (designed for symmetric positive definite systems)
%
%   'CGIters' - maximum number of conjugate gradient iterations allowed {5}
%
%   'CGTolType' - CG stopping tolerance type used ('quadratic')
%     'quadratic' :    || R || / || G || <  min(0.5,|| G ||)
%     'superlinear' :  || R || / || G || <  min(0.5,sqrt(|| G ||))
%     'fixed' :        || R || < CGTol
%   where R is the residual and G is the gradient of FUN at X.
%
%   'CGTol' - CG stopping tolerance when CGTolType is 'fixed' {1e-6}
%
%   'HessVecFDStep' - Hessian vector product finite difference step {1e-10}
%     0 - use the default step given in HESSVEC_FD
%     >0 - fixed value to use at the difference step 
%
%   PARAMS = TN('defaults') returns a structure containing the 
%   default parameters for the particular Poblano method. 
%
%
%   Examples 
%  
%   Suppose the function and gradient of the objective function are
%   specified in an mfile named mysin.m:
%
%     function [f,g]=example1(x,a)
%     if nargin < 2, a = 1; end
%     f = sin(a*x);
%     g = a*cos(a*x);
%
%   We can call the optimization method (using its default
%   parameters) using the command:
%
%     out = tn(@(x) example1(x,3), pi/4);
%
%   To change a parameter, we can specify a param/value input pair
%   as follows:
%
%     out = tn(@(x) example1(x,3), pi/4, 'Display', 'final');
%
%   Alternatively, we can use a structure to define the parameters:
%  
%     params.MaxIters = 2;
%     out = tn(@(x) example1(x,3), pi/4, params);
%
%   See also POBLANO_OUT, POBLANO_PARAMS, HESSVEC_FD, FUNCTION_HANDLE.
%
%MATLAB Poblano Toolbox.
%Copyright 2009-2012, Sandia Corporation.

%% Parse parameters

% Create parser
params = inputParser;

% Set Poblano parameters
params = poblano_params(params);

% Set parameters for this method
params.addParamValue('CGIters',5,@(x) x > 0);
params.addParamValue('CGTolType','quadratic',@(x) ismember(x,{'quadratic','superlinear','fixed'}));
params.addParamValue('CGTol',1e-6, @(x) x > 0);
params.addParamValue('HessVecFDStep',1e-10, @(x) x >= 0);
params.addParamValue('CGSolver','symmlq',@(x) ismember(x,{'symmlq','pcg'}));

% Parse input
params.parse(varargin{:});

%% Check input arguments
if (nargin == 1) && isequal(FUN,'defaults') && (nargout == 1)
    out = params.Results;
    return;
elseif (nargin < 2)
    error('Error: invalid input arguments');
end

%% Initialize
xk = x0;
[fk,gk] = feval(FUN,xk);
out = poblano_out(xk,fk,gk,1,params);
ak = 1.0;

cgIters = params.Results.CGIters;
tolType = params.Results.CGTolType;
cgTol = params.Results.CGTol;
hessvecFDstep = params.Results.HessVecFDStep;
solver = params.Results.CGSolver;

%% Main loop
while out.ExitFlag == -1

    % Compute step direction, pk
    ngk = norm(gk);
    switch tolType
        case 'quadratic'
            cgTol = min(1-eps,min(0.5,ngk)*ngk);       % quadratic convergence
        case 'superlinear'
            cgTol = min(1-eps,min(0.5,sqrt(ngk))*ngk); % superlinear convergence
        otherwise
    end
    
    % Setup Hessian vector product finite difference approximation
    if hessvecFDstep > 0
        hv = @(v) hessvec_fd(v,FUN,xk,gk,hessvecFDstep);
    else
        hv = @(v) hessvec_fd(v,FUN,xk,gk);
    end
    
    % Keep track of number of function calls are made using hessvec_fd
    global nfev_hessvec_fd;
    nfev_hessvec_fd = 0;

    % Compute step direction, pk
    [pk,cg_flag,cg_relres,ncgfev]= feval(solver, hv, -gk, cgTol, cgIters);
    
    % Compute step length
    [xk,fk,gk,ak,lsinfo,nfev] = poblano_linesearch(FUN,xk,fk,gk,ak,pk,params.Results);
    if (lsinfo ~= 1) 
        if strcmp(params.Results.Display, 'iter')
            fprintf(1,[mfilename,': line search warning = %d\n'],lsinfo);
        end
        pk = -gk;
        [xk,fk,gk,ak,lsinfo,nfev] = poblano_linesearch(FUN,xk,fk,gk,ak,pk,params.Results);
    end
    
    % Update counts, check exit conditions, etc.
    out = poblano_out(xk,fk,gk,nfev+nfev_hessvec_fd,params,out);
end

clear nfev_hessvec_fd;
