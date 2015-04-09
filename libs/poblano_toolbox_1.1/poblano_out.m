function out = poblano_out(x,f,g,nfev,params,out)
%POBLANO_OUT   Standard output parameters for the Poblano Toolbox.
%
%   Standard ouptut parameters for Poblano are as follows:
%
%     OUT.X        : solution (i.e., best point found so far)
%     OUT.F        : function value at solution
%     OUT.G        : gradient at solution
%     OUT.Params   : input parameters
%     OUT.FuncEvals: number of function evaluations 
%     OUT.Iters    : number of iterations 
%     OUT.ExitFlag : termination flag
%      -1 = no termination condition has been statisfied
%       0 = successful termination (based on StopTol)
%       1 = maximum number of iterations exceeded
%       2 = maximum number of function values exceeded
%       3 = relative change in function value < RelFuncTol
%       4 = NaNs found in f, g, or norm(g)
%     OUT.ExitDescription : text description of termination flag
%
%   Additional fields may be included if a trace is requested.
%
%   This method should not be called directly, but only by Poblano
%   optimization algorithms.
%
%   See also POBLANO_PARAMS.
%
%MATLAB Poblano Toolbox.
%Copyright 2009-2012, Sandia Corporation.

%% Check if this is the first call
if ~exist('out','var')
    out.Params = params;
    out.ExitFlag = -1;
    out.ExitDescription = '';
    out.X = x;
    out.F = f;
    out.G = g;
    out.FuncEvals = nfev;
    out.Iters = 0;
    if params.Results.TraceX, out.TraceX=[]; end
    if params.Results.TraceFunc, out.TraceFunc=[]; end
    if params.Results.TraceRelFunc, out.TraceRelFunc=[]; end
    if params.Results.TraceGrad, out.TraceGrad=[]; end
    if params.Results.TraceGradNorm, out.TraceGradNorm=[]; end
    if params.Results.TraceFuncEvals, out.TraceFuncEvals=[]; end
else
    oldf = out.F;
    if (f <= out.F)
        out.X = x;
        out.F = f;
        out.G = g;
    end
    if abs(oldf) < eps
        relfit = abs(f - oldf);
    else
        relfit = abs((f - oldf) / oldf);
    end
    out.FuncEvals = out.FuncEvals + nfev;  
    out.Iters = out.Iters + 1;
end

%% Initialize variables
nx = length(x);
g2norm = norm(g);
g2normnx = g2norm/nx;

%% Record trace info
if params.Results.TraceX
    out.TraceX(:,end+1) = x; 
end

if params.Results.TraceFunc 
    out.TraceFunc(end+1) = f; 
end

if  params.Results.TraceFunc && params.Results.TraceRelFunc 
    if (out.Iters > 0)
        out.TraceRelFunc(end+1) = relfit;
    end
end

if params.Results.TraceGrad 
    out.TraceGrad(:,end+1) = g; 
end

if params.Results.TraceGradNorm
    out.TraceGradNorm(end+1) = g2norm; 
end

if params.Results.TraceFuncEvals
    out.TraceFuncEvals(end+1) = nfev;
end

%% Check termination conditions
if g2normnx < params.Results.StopTol 
    % solution found
    out.ExitFlag = 0;
    out.ExitDescription = 'Successful termination based on StopTol';
elseif out.Iters >= params.Results.MaxIters
    % maximum iterations exceeded
    out.ExitFlag = 1;
    out.ExitDescription = 'Maximum number of iterations exceeded';
elseif out.FuncEvals >= params.Results.MaxFuncEvals, 
    % maximum function evaluations exceeded
    out.ExitFlag = 2; 
    out.ExitDescription = 'Maximum number of function evaluations exceeded';
elseif (out.Iters > 0) && (relfit <= params.Results.RelFuncTol)
    % relative fit tolerance is reached
    out.ExitFlag=3;            
    out.ExitDescription = 'Relative change in F < RelFuncTol';
elseif isnan(f) || sum(isnan(g)) || isnan(g2norm)
    % new point found results in NaNs
    out.ExitFlag = 4;
    out.ExitDescription = 'NaNs found in F, G, or norm(G)';
end

%% Display iteration information

% First Iteration
if (out.Iters == 0) && ~strcmp(params.Results.Display, 'off')
    fprintf(1,' Iter  FuncEvals       F(X)          ||G(X)||/N        \n');
    fprintf(1,'------ --------- ---------------- ----------------\n');
end
    
% Iteration info
if strcmp(params.Results.Display, 'iter') && (mod(out.Iters,params.Results.DisplayIters)==0 || out.ExitFlag>=0)
    fprintf(1,'%6d %9d %16.8f %16.8f\n', out.Iters, out.FuncEvals, ...
            f, g2normnx);
end

% Final Iteration
if ((out.ExitFlag >= 0) && strcmp(params.Results.Display, 'final'))
%    if strcmp(params.Results.Display, 'iter')
%        fprintf(1,'------ --------- ---------------- ----------------\n');
%    end
    fprintf(1,'%6d %9d %16.8f %16.8f\n', out.Iters, out.FuncEvals, ...
            out.F, g2normnx);
end
