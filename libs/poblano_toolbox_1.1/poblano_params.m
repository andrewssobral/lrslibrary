function params = poblano_params(params)
%POBLANO_PARAMS Standard input parameters for the Poblano Toolbox.
%
%   Standard input parameters for Poblano are as follows, with defaults in
%   parentheses:
%
%   Display - controls printed output {'iter'}
%     'iter' : display every iteration
%     'final': display only after final iteration
%     'off'  : no display
%
%   DisplayIters - number of iterations to display printed output {1}
%
%   MaxIters - maximum number of iterations {1000}
%
%   MaxFuncEvals - maximum number of function evaluations {10000}
%
%   StopTol - gradient norm stopping tolerance {1e-5}, i.e., the
%   method stops when the norm of the gradient is less than StopTol
%   times the number of variables.
%
%   RelFuncTol - relative function value change stopping tolerance {1e-5}, 
%   i.e., the method stops when the relative change of the function value 
%   from one iteration to the next is less than RelFuncTol.
%
%   TraceX, TraceFunc, TraceRelFunc, TraceGrad, TraceGradNorm,
%   TraceFuncEvals {false} -  booleans specifying whether or not to
%   record X, F(X), G(X), ||G(X)||, and number of function
%   evaluations, respectively, at each iteration. 
%
%   LineSearch_method {'more-thuente'} - string specifying the line search
%   method to use
%     'more-thuente': find a step which satisfies a sufficient decrease
%                     condition and a curvature condition (CVSRCH)
%
%   See POBLANO_LINESEARCH for details on the line search parameters for 
%   each available method. The available parameters:
%     LineSearch_xtol {1e-15}
%     LineSearch_ftol {1e-4}
%     LineSearch_gtol {1e-2}
%     LineSearch_stpmin {1e-15}
%     LineSearch_stpmax {1e15}
%     LineSearch_maxfev {20}
%     LineSearch_initialstep {1}
%
%   This function should not be called directly, but only by Poblano
%   optimization algorithms.
%
%
%   Examples
%
%   The default options can be retrieved from any Poblano method by
%
%   %  params = POBLANO_METHOD('defaults')
%      params = ncg('defaults')
%
%   Also, it is possible to pass in parameters to any Poblano method as
%   a strucure.
%
%   See also POBLANO_LINESEARCH.
%
%MATLAB Poblano Toolbox.
%Copyright 2009-2012, Sandia Corporation.

%% Standard parameters across Poblano methods
params.addParamValue('Display','iter', @(x) ismember(x,{'iter','final','off'}));
params.addParamValue('DisplayIters', 1, @(x) x >= 1);
params.addParamValue('MaxIters', 1000, @(x) x >= 0);
params.addParamValue('MaxFuncEvals', 10000, @(x) x >= 0); 
params.addParamValue('StopTol', 1e-5, @(x) x >= 0);
params.addParamValue('RelFuncTol', 1e-6, @isnumeric);
params.addParamValue('TraceX', false, @islogical);
params.addParamValue('TraceFunc', false, @islogical);
params.addParamValue('TraceRelFunc', false, @islogical);
params.addParamValue('TraceGrad', false, @islogical);
params.addParamValue('TraceGradNorm', false, @islogical);
params.addParamValue('TraceFuncEvals', false, @islogical);
params.addParamValue('LineSearch_method', 'more-thuente', @(x) ismember(x,{'more-thuente'}));
params.addParamValue('LineSearch_initialstep', 1, @(x) x >= 0);
params.addParamValue('LineSearch_xtol', 1e-15, @(x) x >= 0);
params.addParamValue('LineSearch_ftol', 1e-4, @(x) x >= 0);
params.addParamValue('LineSearch_gtol', 1e-2, @(x) x >= 0);
params.addParamValue('LineSearch_stpmin', 1e-15, @(x) x >= 0);
params.addParamValue('LineSearch_stpmax', 1e15, @(x) x >= 0);
params.addParamValue('LineSearch_maxfev', 20, @(x) x >= 0);
