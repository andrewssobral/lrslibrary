function [alpha,output] = ls_mt(f,g,z,p,state,options)
%LS_MT Strong Wolfe line search by More-Thuente.
%   [alpha,output] = ls_mt(f,g,z,p) attempts to find a point along the line
%   z+alpha*p that satisfies the strong Wolfe conditions for the function f
%   and its gradient g. In other words, a step alpha is computed so that
%   both
%
%      f(z+alpha*p) <= f(z)+c1*alpha*real(g(z)'*p)           (suff. decr.)
%      abs(real(g(z+alpha*p)'*p)) <= c2*abs(real(g(z)'*p))   (curv. cond.)
%
%   where p is a descent direction and 0 < c1 < c2 < 1, are satisfied. If g
%   is the empty matrix [], it is approximated with finite differences. The
%   variables z and p may be a scalar, vector, matrix, tensor or even a
%   cell array of tensors and their contents may be real or complex.
%
%   If f(x) is a function of real variables x, the function g(x) should
%   compute the partial derivatives of f with respect to the real variables
%   x, i.e. g(xk) := df(xk)/dx. If f(z) is a function of complex variables
%   z, the function g(z) should compute two times the partial derivative
%   of f with respect to conj(z) (treating z as constant), i.e. g(z) :=
%   2*df(z)/d(conj(z)) = 2*conj(df(z)/dz). The output of the function g(z)
%   may have the same structure as z (although this is not necessary).
%
%   [alpha,output] = ls_mt(f,g,z,p,state) optionally supplies the function
%   value and gradient at z in state. The fields state.fval(end) and
%   state.grad should equal f(z) and g(z), respectively. If state is the
%   empty matrix [], or if state.fval or state.grad are missing, these
%   fields are evaluated using f and g. The field state.grad may have the
%   same structure as z, but may also be given in vectorized form.
%
%   The structure output returns additional information:
%      
%      output.dgrad   - The directional derivative real(output.grad'*p).
%      output.fevals  - The amount of function/gradient calls.
%      output.fval    - The value f(z+alpha*p).
%      output.grad    - The vectorized output of g(z+alpha*p).
%      output.info    - The circumstances under which the procedure
%                       terminated:
%                          1: Strong Wolfe conditions are satisfied.
%                          2: Maximum number of function calls reached.
%                          3: Step is equal to the minimum step size.
%                          4: Step is equal to the maximum step size.
%                          5: The minimum has been bracketed and the
%                             bracket width is less than the specified
%                             tolerance.
%                          6: Rounding errors prevent progress.
%      output.message - The termination message.
%      output.z       - The new iterate z+alpha*p.
%
%   ls_mt(f,g,z,p,state,options) may be used to set the following
%   options:
%
%      options.alpha = 1        - The initial step size.
%      options.c1 = 1e-4        - The sufficient decrease parameter.
%      options.c2 = 0.9         - The curvature condition parameter.
%      options.MaxFunEvals = 30 - The maximum amount of function calls.
%      options.MaxStep = inf    - The maximum step size.
%      options.MinStep = 0      - The minimum step size.
%      options.TolX = 1e-9      - The relative tolerance for the bracket
%                                 width.

%   Authors: Original MINPACK-2 implementation by Brett M. Averick,
%            Richard G. Carter and Jorge J. More.
%            Translated to MATLAB and adapted for functions of complex
%            variables by Laurent Sorber (Laurent.Sorber@cs.kuleuven.be).
%
%   References:
%   [1] J. J. More, D. J. Thuente, "Line search algorithms with guaranteed
%       sufficient decrease", ACM Transactions on Mathematical Software,
%       Vol. 20, No. 3, 1994, pp. 286-307.
%   [2] L. Sorber, M. Van Barel, L. De Lathauwer, "Unconstrained
%       optimization of real functions in complex variables", SIAM J. Opt.,
%       Vol. 22, No. 3, 2012, pp. 879-898.

% Check the options structure.
if nargin < 6, options = struct; end
if ~isfield(options,'alpha'), options.alpha = 1; end
if ~isfield(options,'c1'), options.c1 = 1e-4; end
if ~isfield(options,'c2'), options.c2 = 0.9; end
if ~isfield(options,'MaxFunEvals'), options.MaxFunEvals = 30; end
if ~isfield(options,'MaxStep'), options.MaxStep = inf; end
if ~isfield(options,'MinStep'), options.MinStep = 0; end
if ~isfield(options,'TolX'), options.TolX = 1e-9; end

% Initialize output.
output.fevals = 0;
output.info = NaN;

% Convert input to MINPACK-2 format.
if isfield(state,'fval') && ~isempty(state.fval)
    fk = state.fval(end);
else
    fk = f(z);
end
if isfield(state,'grad') && ~isempty(state.grad)
    gk = state.grad;
else
    if ~isa(g,'function_handle') && isempty(g)
        gk = deriv(f,z,fk);
    else
        gk = g(z);
    end
end
dim = structure(z);
z = serialize(z);
p = serialize(p);
gk = serialize(gk);
stp = options.alpha;
% Initially: f = fk. In the loop: f = fstp.
dgk = real(gk'*p); % Initially: g = dgk. In the loop: g = dgstp.
ftol = options.c1;
gtol = options.c2;
xtol = options.TolX;
stpmin = options.MinStep;
stpmax = options.MaxStep;

% subroutine dcsrch(stp,f,g,ftol,gtol,xtol,task,stpmin,stpmax,isave,dsave)
p5 = 0.5;
p66 = 0.66;
xtrapl = 1.1;
xtrapu = 4.0;

% Check the input arguments for errors.
if stp < stpmin
    error('ls_mt:alpha', ...
          'The initial step size is smaller than the minimum step size.');
end
if stp > stpmax
    error('ls_mt:alpha', ...
          'The initial step size is larger than the maximum step size.');
end
if dgk > 0
    error('ls_mt:gk', ...
          'Gradient at z is not a descent direction.');
end
if ftol < 0
    error('ls_mt:c1', ...
          'The parameter c1 should be greater than zero.');
end
if gtol < 0
    error('ls_mt:c2', ...
          'The parameter c2 should be greater than zero.');
end
if xtol < 0
    error('ls_mt:TolX','TolX should be greater than zero.');
end
if stpmin < 0
    error('ls_mt:MinStep', ...
          'The minimum step size should be greater than zero.');
end
if stpmax < stpmin
    error('ls_mt:MaxStep', ...
         ['The maximum step size should be greater than the minimum ' ...
          'step size.']);
end

% Initialize local variables.
brackt = false;
stage = 1;
finit = fk;
dginit = dgk;
dgtest = ftol*dginit;
width = stpmax-stpmin;
width1 = width/p5;

% The variables stx, fx, dgx contain the values of the step, function, and
% derivative at the best step. The variables sty, fy, dgy contain the value
% of the step, function, and directional derivative at sty. The variables
% stp, fk, dgk contain the values of the step, function, and directional
% derivative at stp.
stx = 0;
fx = finit;
dgx = dginit;
sty = 0;
fy = finit;
dgy = dginit;
stmin = 0;
stmax = stp+xtrapu*stp;

while true

    % Obtain another function and derivative.
    zstp = deserialize(z+stp*p,dim);
    fstp = f(zstp);
    if ~isa(g,'function_handle') && isempty(g)
        gstp = serialize(deriv(f,zstp,fstp));
    else
        gstp = serialize(g(zstp));
    end
    dgstp = real(gstp'*p);
    output.fevals = output.fevals+1;

    % If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the algorithm
    % enters the second stage.
    ftest = finit+stp*dgtest;
    if stage == 1 && fstp <= ftest && dgstp >= 0
        stage = 2;
    end

    % Test for warnings and convergence.
    if brackt && (stp <= stmin || stp >= stmax)
        output.info = 6;
        output.message = 'Rounding errors prevent progress.';
    end
    if brackt && stmax-stmin <= xtol*stmax
        output.info = 5;
        output.message = ['The minimum has been bracketed and the ' ...
                          'bracket width is less than the specified ' ...
                          'tolerance.'];
    end
    if stp == stpmax && fstp <= ftest && dgstp <= dgtest
        output.info = 4;
        output.message = 'Step is equal to the maximum step size.';
    end
    if stp == stpmin && (fstp > ftest || dgstp >= dgtest)
        output.info = 3;
        output.message = 'Step is equal to the minimum step size.';
    end
    if output.fevals >= options.MaxFunEvals
        output.info = 2;
        output.message = 'Maximum number of function evaluations reached.';
    end
    if fstp <= ftest && abs(dgstp) <= gtol*(-dginit)
        output.info = 1;
        output.message = 'Strong Wolfe conditions are satisfied.';
    end

    % Test for termination. 
    if ~isnan(output.info)
        alpha = stp;
        output.z = zstp;
        output.fval = fstp;
        output.grad = gstp;
        output.dgrad = dgstp;
        return;
    end

    % A modified function is used to predict the step during the first
    % stage if a lower function value has been obtained but the decrease is
    % not sufficient.
    if stage == 1 && fstp <= fx && fstp > ftest

        % Define the modified function and derivative values.
        fstpm = fstp-stp*dgtest;
        fxm = fx-stx*dgtest;
        fym = fy-sty*dgtest;
        dgstpm = dgstp-dgtest;
        dgxm = dgx-dgtest;
        dgym = dgy-dgtest;

        % Call dcstep to update stx, sty, and to compute the new step.
        %call dcstep(stx,fxm,gxm,sty,fym,gym,stp,fm,gm,brackt,stmin,stmax)
        [fxm,dgxm,fym,dgym] = dcstep(fxm,dgxm,fym,dgym,fstpm,dgstpm, ...
            stmin,stmax);

        % Reset the function and derivative values for fstp.
        fx = fxm+stx*dgtest;
        fy = fym+sty*dgtest;
        dgx = dgxm+dgtest;
        dgy = dgym+dgtest;

    else

        % Call dcstep to update stx, sty, and to compute the new step.
        %call dcstep(stx,fx,gx,sty,fy,gy,stp,f,g,brackt,stmin,stmax)
        [fx,dgx,fy,dgy] = dcstep(fx,dgx,fy,dgy,fstp,dgstp,stmin,stmax);

    end

    % Decide if a bisection step is needed.
    if brackt == true
        if abs(sty-stx) >= p66*width1
            stp = stx+p5*(sty-stx);
        end
        width1 = width;
        width = abs(sty-stx);
    end

    % Set the minimum and maximum steps allowed for stp. 
    if brackt == true
        stmin = min(stx,sty);
        stmax = max(stx,sty);
    else
        stmin = stp+xtrapl*(stp-stx);
        stmax = stp+xtrapu*(stp-stx);
    end

    % Force the step to be within the bounds stpmax and stpmin.
    stp = max(stp,stpmin);
    stp = min(stp,stpmax);

    % If further progress is not possible, let stp be the best point
    % obtained during the search.
    if (brackt && (stp <= stmin || stp >= stmax) || ...
       (brackt && stmax-stmin <= xtol*stmax))
        stp = stx;
    end

end

function [fx,dx,fy,dy] = dcstep(fx,dx,fy,dy,fp,dp,stpmin,stpmax)
% subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,stpmin,stpmax)

sgnd = dp*(dx/abs(dx));

% First case: A higher function value. The minimum is bracketed.
% If the cubic step is closer to stx than the quadratic step, the cubic
% step is taken, otherwise the average of the cubic and quadratic steps
% is taken.
if fp > fx
    theta = 3*(fx-fp)/(stp-stx)+dx+dp;
    s = max([abs(theta) abs(dx) abs(dp)]);
    gamma = s*sqrt((theta/s)^2-(dx/s)*(dp/s));
    if stp < stx, gamma = -gamma; end
    o = (gamma-dx)+theta;
    q = ((gamma-dx)+gamma)+dp;
    r = o/q;
    stpc = stx+r*(stp-stx);
    stpq = stx+((dx/((fx-fp)/(stp-stx)+dx))/2)*(stp-stx);
    if abs(stpc-stx) < abs(stpq-stx)
        stpf = stpc;
    else
        stpf = stpc+(stpq-stpc)/2;
    end
    brackt = true;

% Second case: A lower function value and derivatives of opposite sign.
% The minimum is bracketed. If the cubic step is farther from stp than 
% the secant step, the cubic step is taken, otherwise the secant step 
% is taken.
elseif sgnd < 0
    theta = 3*(fx-fp)/(stp-stx)+dx+dp;
    s = max([abs(theta) abs(dx) abs(dp)]);
    gamma = s*sqrt((theta/s)^2-(dx/s)*(dp/s));
    if stp > stx, gamma = -gamma; end
    o = (gamma-dp)+theta;
    q = ((gamma-dp)+gamma)+dx;
    r = o/q;
    stpc = stp+r*(stx-stp);
    stpq = stp+(dp/(dp-dx))*(stx-stp);
    if abs(stpc-stp) > abs(stpq-stp)
        stpf = stpc;
    else
        stpf = stpq;
    end
    brackt = true;

% Third case: A lower function value, derivatives of the same sign, and
% the magnitude of the derivative decreases.
elseif abs(dp) < abs(dx)
    % The cubic step is computed only if the cubic tends to infinity in
    % the direction of the step or if the minimum of the cubic is
    % beyond stp. Otherwise the cubic step is defined to be the secant
    % step.
    theta = 3*(fx-fp)/(stp-stx)+dx+dp;
    s = max([abs(theta) abs(dx) abs(dp)]);
    % The case gamma = 0 only arises if the cubic does not tend to 
    % infinity in the direction of the step.
    gamma = s*sqrt(max(0,(theta/s)^2-(dx/s)*(dp/s)));
    if (stp > stx), gamma = -gamma; end
    o = (gamma-dp)+theta;
    q = (gamma+(dx-dp))+gamma;
    r = o/q;
    if r < 0 && gamma ~= 0
        stpc = stp+r*(stx-stp);
    elseif stp > stx
        stpc = stpmax;
    else
        stpc = stpmin;
    end
    stpq = stp+(dp/(dp-dx))*(stx-stp);
    if brackt == true
        % A minimizer has been bracketed. If the cubic step is closer
        % to stp than the secant step, the cubic step is taken, 
        % otherwise the secant step is taken.
        if abs(stpc-stp) < abs(stpq-stp)
            stpf = stpc;
        else
            stpf = stpq;
        end
        if stp > stx
            stpf = min(stp+p66*(sty-stp),stpf);
        else
            stpf = max(stp+p66*(sty-stp),stpf);
        end
    else
        % A minimizer has not been bracketed. If the cubic step is 
        % farther from stp than the secant step, the cubic step is
        % taken, otherwise the secant step is taken.
        if abs(stpc-stp) > abs(stpq-stp)
            stpf = stpc;
        else
            stpf = stpq;
        end
        stpf = min(stpmax,stpf);
        stpf = max(stpmin,stpf);
    end

% Fourth case: A lower function value, derivatives of the same sign,
% and the magnitude of the derivative does not decrease. If the minimum
% is not bracketed, the step is either stpmin or stpmax, otherwise the
% cubic step is taken.
else
    if brackt == true
        theta = 3*(fp-fy)/(sty-stp)+dy+dp;
        s = max([abs(theta) abs(dy) abs(dp)]);
        gamma = s*sqrt((theta/s)^2-(dy/s)*(dp/s));
        if stp > sty, gamma = -gamma; end
        o = (gamma-dp)+theta;
        q = ((gamma-dp)+gamma)+dy;
        r = o/q;
        stpc = stp+r*(sty-stp);
        stpf = stpc;
    elseif stp > stx
        stpf = stpmax;
    else
        stpf = stpmin;
    end

end

% Update the interval which contains a minimizer.
if fp > fx
    sty = stp;
    fy = fp;
    dy = dp;
else
    if sgnd < 0
        sty = stx;
        fy = fx;
        dy = dx;
    end
    stx = stp;
    fx = fp;
    dx = dp;
end

% Compute the new step.
stp = stpf;

end

end

function z = deserialize(z,dim)
    if iscell(dim)
        v = z; z = cell(size(dim));
        s = cellfun(@(s)prod(s(:)),dim(:)); o = [0; cumsum(s)];
        for i = 1:length(s), z{i} = reshape(v(o(i)+(1:s(i))),dim{i}); end
    elseif ~isempty(dim)
        z = reshape(z,dim);
    end
end

function z = serialize(z)
    if iscell(z)
        s = cellfun(@numel,z(:)); o = [0; cumsum(s)];
        c = z; z = zeros(o(end),1);
        for i = 1:length(s), ci = c{i}; z(o(i)+(1:s(i))) = ci(:); end
    else
        z = z(:);
    end
end

function dim = structure(z)
    if iscell(z)
        dim = cellfun(@size,z,'UniformOutput',false);
    else
        dim = size(z);
        if numel(z) == dim(1), dim = []; end
    end
end
