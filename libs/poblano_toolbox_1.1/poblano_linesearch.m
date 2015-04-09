function [x,f,g,a,lsinfo,lsnfev] = poblano_linesearch(FUN,x0,f0,g0,a0,d0,params)
%POBLANO_LINESEARCH Line search methods in the Poblano Toolbox.
%
%  [X,F,G,A] = POBLANO_LINESEARCH(FUN,X0,F0,G0,A0,D0,PARAMS) performs a
%  line search to find a point X such that FUN(X0+A*D0) is minimized. The
%  method returns X = X0+A*D0 with the function (F) and gradient (G) at
%  that point, along with the step length (A) found by the line search. The
%  parameters used in the line search are specified in PARAMS (see
%  POBLANO_PARAMS for more details).
% 
%  If PARAMS.LineSearch_initialstep equals 0, the initial step used in the
%  line search is A0, which is 1 for the first iteration and the step
%  length computed in the previous iteration for all subsequent iterations
%  If PARAMS.LineSearch_initialstep > 0, the initial step is the value of
%  PARAMS.LineSearch_initialstep.
%
%  [X,F,G,A,LSINFO] = POBLANO_LINESEARCH(FUN,...) also provides an exit
%  code from the line search method. See the particular line search method
%  for an explanation of the exit codes.
%
%  [X,F,G,A,LSINFO,LSNFEV] = POBLANO_LINESEARCH(FUN,...) also provides the
%  number of function/gradient calculations performed during the line search.
%
%  Parameters (PARAMS)
%  
%  PARAMS.LineSearch_method
%
%    'more-thuente' : More-Thuente line search from MINPACK, adapted for 
%                     Matlab by Dianne O'Leary
%
%      See CVSRCH for details on the use of parameters. Below are how the
%      Poblano parameters map to the CVSRCH parameters:
%
%      LineSearch_xtol: xtol
%      LineSearch_ftol: ftol
%      LineSearch_gtol: gtol
%      LineSearch_stpmin: stpmin
%      LineSearch_stpmax: stpmax
%      LineSearch_maxfev: maxfev
%
%  PARAMS.LineSearch_initialstep
%    0  : use step length provided in A0
%    >0 : use step length provided in PARAMS.Linesearch_initialstep 
%
%   See also POBLANO_PARAMS, CVSRCH.
%
%MATLAB Poblano Toolbox.
%Copyright 2009-2012, Sandia Corporation.

if nargin < 7
    error('POBLANO_LINESEARCH: too few arguments');
end

switch (params.LineSearch_method)
    case 'more-thuente'
        minfun = 'cvsrch';
end

% Check whether user specified an initial step
if params.LineSearch_initialstep > 0
    a0 = params.LineSearch_initialstep;
end

[x,f,g,a,lsinfo,lsnfev] = feval(minfun,FUN,x0,f0,g0,a0,d0,params);
