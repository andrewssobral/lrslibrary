function [z,output] = minf_lbfgsdl(f,g,z0,options)
%MINF_LBFGSDL Minimize a function by L-BFGS with dogleg trust region.
%   [z,output] = minf_lbfgsdl(f,g,z0) starts at z0 and attempts to find a
%   local minimizer of the real-valued function f(z). The input variables z
%   may be a scalar, vector, matrix, tensor or even a (nested) cell array
%   of tensors and its contents may be real or complex.
%
%   If f(x) is a function of real variables x, the function g(x) should
%   compute the partial derivatives of f with respect to the real variables
%   x, i.e. g(xk) := df(xk)/dx. If f(z) is a function of complex variables
%   z, the function g(z) should compute two times the partial derivative
%   of f with respect to conj(z) (treating z as constant), i.e. g(zk) :=
%   2*df(zk)/d(conj(z)) = 2*conj(df(zk)/dz). If g is the empty matrix [],
%   the real gradient or scaled conjugate cogradient is approximated with
%   finite differences. The output of the function g(z) may have the same
%   structure as z (although this is not necessary). The structure output
%   returns additional information:
%
%      output.alpha      - The plane search step lengths in every
%                          iteration, if a plane search is selected.
%      output.delta      - The trust region radius at every step attempt.
%      output.fevals     - The total number of function calls.
%      output.fval       - The value of the objective function f in every
%                          iteration.
%      output.gevals     - The total number of gradient calls.
%      output.info       - The circumstances under which the procedure
%                          terminated:
%                             1: Objective function tolerance reached.
%                             2: Step size tolerance reached.
%                             3: Maximum number of iterations reached.
%      output.infops     - The circumstances under which the plane search
%                          terminated in every iteration.
%      output.iterations - The number of iterations.
%      output.relfval    - The difference in objective function value
%                          between every two successive iterates, relative
%                          to its initial value.
%      output.relstep    - The step size relative to the norm of the 
%                          current iterate in every iteration.
%      output.rho        - The trustworthiness at every step attempt.
%
%   minf_lbfgsdl(f,g,z0,options) may be used to set the following options:
%
%      options.Delta =            - The initial trust region radius.
%      0.3*max(1,norm(z0))
%      options.Display = 10       - Displays the objective function value,
%                                   its difference with the previous
%                                   iterate relative to the first iterate
%                                   and the relative step size each
%                                   options.Display iterations. Set to 0 to
%                                   disable.
%      options.M =                - The number of updates to store.
%      min(30,length(z0))
%      options.MaxIter = 500      - The maximum number of iterations.
%      options.PlaneSearch        - The plane search used to minimize the
%      = false                      objective function in the plane spanned
%                                   by the steepest descent direction and
%                                   the Gauss-Newton step. Disables dogleg
%                                   trust region strategy. The method
%                                   should have the function signature
%                                   options.PlaneSearch(F,dF,z,p1,p2, ...
%                                   state,options.PlaneSearchOptions).
%      options.PlaneSearchOptions - The options structure passed to the
%                                   plane search search routine.
%      options.TolFun = 1e-6      - The tolerance for output.relfval.
%      options.TolX = 1e-8        - The tolerance for output.relstep.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Unconstrained
%       optimization of real functions in complex variables", SIAM J. Opt.,
%       Vol. 22, No. 3, 2012, pp. 879-898.

% Store the structure of the input space and evaluate the gradient.
dim = structure(z0);
fval = f(z0);
if ~isa(g,'function_handle') && isempty(g)
    grad = serialize(deriv(f,z0,fval));
else
    grad = serialize(g(z0));
end
z0 = serialize(z0);

% Check the options structure.
if nargin < 4, options = struct; end
if ~isfield(options,'Delta'), options.Delta = 0.3*max(1,norm(z0)); end
if ~isfield(options,'Display'), options.Display = 10; end
if ~isfield(options,'M'), options.M = min(30,length(z0)); end
if ~isfield(options,'MaxIter'), options.MaxIter = 500; end
if ~isfield(options,'PlaneSearch'), options.PlaneSearch = false; end
if ~isfield(options,'PlaneSearchOptions')
    options.PlaneSearchOptions = struct;
end
if ~isfield(options,'TolFun'), options.TolFun = 1e-6; end
if ~isfield(options,'TolX'), options.TolX = 1e-8; end

% Initialize the algorithm.
S = zeros(numel(z0),options.M);
Y = zeros(numel(z0),options.M);
a = zeros(1,options.M);
r = zeros(1,options.M);
m = 0;
midx = 0;

% L-BFGS with dogleg trust region.
output.alpha = [];
output.delta = options.Delta;
output.fevals = 1;
output.fval = fval;
output.gevals = 1;
output.info = false;
output.infops = [];
output.iterations = 0;
output.relfval = [];
output.relstep = [];
output.rho = [];
while ~output.info

    % Compute the quasi-Newton step pqn = -H*grad.
    pqn = -grad;
    for i = 1:m
        a(i) = r(midx(i))*real(S(:,midx(i))'*pqn);
        pqn = pqn-a(i)*Y(:,midx(i));
    end
    if m > 0
        y1 = Y(:,midx(1));
        y1y1 = y1'*y1;
        gamma = y1y1*r(midx(1));
        pqn = 1/gamma*pqn;
    end
    for i = m:-1:1
        b = r(midx(i))*real(Y(:,midx(i))'*pqn);
        pqn = pqn+(a(i)-b)*S(:,midx(i));
    end
    
    % Approximate the Cauchy point pcp = -alpha*grad, where alpha is equal
    % to arg min m(-alpha*grad) and m(p) is a second order model of f at z.
    gg = grad'*grad;
    if m == 0
        alpha = 1;
    else
        s1 = S(:,midx(1));
        gBg = gg-real(grad'*s1)^2/(s1'*s1)+real(grad'*y1)^2/y1y1;
        gBg = gamma*gBg;
        alpha = gg/gBg;
    end
    
    % Plane search in the plane spanned by {pqn,pcp}.
    if isa(options.PlaneSearch,'function_handle')
        if output.iterations == 0, pcp = zeros(size(grad));
        else pcp = -alpha*grad; end
        state = output; state.grad = grad;
        [alpha,outputps] = options.PlaneSearch( ...
            f,g,deserialize(z0,dim),deserialize(pqn,dim), ...
            deserialize(pcp,dim),state,options.PlaneSearchOptions);
        output.alpha(:,end+1) = alpha;
        if length(alpha) < 3, alpha(3) = 1; end
        p = alpha(1)*pqn+alpha(2)*pcp;
        z = deserialize(alpha(3)*(z0+p),dim);
        relstep = norm(p)/norm(z0); if isnan(relstep), relstep = 0; end
        if isfield(outputps,'fval'), fval = outputps.fval;
        else fval = f(z); end
        if isfield(outputps,'info')
            output.infops(end+1) = outputps.info;
        end
        rho = 1;
    else
        rho = -inf;
    end

    % Dogleg trust region.
    normpqn = norm(pqn);
    while rho <= 0

        % Compute the dogleg step p.
        delta = output.delta(end);
        if normpqn <= delta
            p = pqn;
            dfval = -0.5*real(grad'*pqn);
        elseif abs(alpha)*sqrt(gg) >= delta
            p = (-delta/sqrt(gg))*grad;
            dfval = delta*(sqrt(gg)-0.5*delta/alpha);
        else
            bma = pqn+alpha*grad; bmabma = bma'*bma;
            a = -alpha*grad; aa = alpha^2*gg;
            c = real(a'*bma);
            if c <= 0
                beta = (-c+sqrt(c^2+bmabma*(delta^2-aa)))/bmabma;
            else
                beta = (delta^2-aa)/(c+sqrt(c^2+bmabma*(delta^2-aa)));
            end
            p = a+beta*bma;
            dfval = 0.5*alpha*(1-beta)^2*gg- ...
                    0.5*beta *(2-beta)*real(grad'*pqn);
        end

        % Compute the trustworthiness rho.
        if dfval > 0
            z = deserialize(z0+p,dim);
            fval = f(z);
            rho = (output.fval(end)-fval)/dfval;
            if isnan(rho), rho = -inf; end
            output.rho(end+1) = rho;
            output.fevals = output.fevals+1;
        end

        % Update trust region radius delta.
        if rho > 0.5
            output.delta(end+1) = max(delta,2*norm(p));
        else
            sigma = (1-0.25)/(1+exp(-14*(rho-0.25)))+0.25;
            if normpqn < sigma*delta && rho < 0
                e = ceil(log2(normpqn/delta)/log2(sigma));
                output.delta(end+1) = sigma^e*delta;
            else
                output.delta(end+1) = sigma*delta;
            end
        end
        
        % Check for convergence.
        relstep = norm(p)/norm(z0); if isnan(relstep), relstep = 0; end
        if rho <= 0 && relstep <= options.TolX
            output.rho(end+1) = rho;
            fval = output.fval(end);
            z = deserialize(z0,dim);
            break;
        end

    end
    
    % Save current state.
    if rho > 0
        z0 = serialize(z);
        grad1 = grad;
    end

    % Evaluate the gradient and update step information.
    if rho > 0
        if isa(options.PlaneSearch,'function_handle') && ...
           output.iterations >= 1 && isfield(outputps,'grad')
            grad = outputps.grad;
        elseif ~isa(g,'function_handle') && isempty(g)
            grad = serialize(deriv(f,z,fval));
        else
            grad = serialize(g(z));
        end
        s = p;
        y = grad-grad1;
        sy = real(y'*s);
        if sy > 0
            m = min(m+1,options.M);
            midx = [midx(1)+1:-1:1,m:-1:midx(1)];
            S(:,midx(1)) = s;
            Y(:,midx(1)) = y;
            r(:,midx(1)) = 1/sy;
        end
    end
    
    % Update the output structure.
    output.fval(end+1) = fval;
    output.gevals = output.gevals+1;
    output.iterations = output.iterations+1;
    output.relfval(end+1) = ...
        abs(diff(output.fval(end:-1:end-1)))/abs(output.fval(1));
    output.relstep(end+1) = relstep;
    if output.relfval(end) <= options.TolFun, output.info = 1; end
    if output.relstep(end) <= options.TolX, output.info = 2; end
    if output.iterations >= options.MaxIter, output.info = 3; end
    
    % Display progress.
    if options.Display > 0 && (output.iterations == 1 || output.info || ...
       mod(output.iterations,options.Display) == 0)
        if output.iterations == 1
            bold = '%s';
            [~,~,~,~,v] = regexp(version('-release'),'([0-9]+)([ab])');
            if usejava('Desktop') && str2double(v{1}{1}) > 2011 || ...
               (str2double(v{1}{1}) == 2011 && strcmpi(v{1}{2},'b'))
                bold = '<strong>%s</strong>';
            end
        end
        if output.iterations == 1 || ...
           mod(output.iterations,15*options.Display) == 0
            fprintf('\n%7s%s','',sprintf(bold,'fval'));
            fprintf('%13s%s','',sprintf(bold,'relfval'));
            fprintf('%10s%s','',sprintf(bold,'relstep'));
            if isa(options.PlaneSearch,'function_handle')
                fprintf('%10s%s','',sprintf(bold,'alpha'));
            else
                fprintf('%10s%s','',sprintf(bold,'delta'));
                fprintf('%8s%s','',sprintf(bold,'rho'));
            end
            fprintf('\n%21s%9s = %4.e %6s = %4.e\n\n','=1/2*norm(F)^2', ...
                    'TolFun',options.TolFun,'TolX',options.TolX);
        end
        if output.iterations == 1
            fprintf('%4i: % 14.8e |\n',0,output.fval(1));
        end
        if isa(options.PlaneSearch,'function_handle')
            stralpha = [repmat('%10.4e ',1,size(output.alpha,1)) '\n'];
            fprintf(['%4i: % 14.8e | %14.8e | %14.8e | ' stralpha], ...
                    output.iterations,output.fval(end), ...
                    output.relfval(end),output.relstep(end), ...
                    abs(output.alpha(:,end)));
        else
            fprintf(['%4i: % 14.8e | %14.8e | %14.8e | ' ...
                     '%10.4e | %10.4e\n'],...
                    output.iterations,output.fval(end), ...
                    output.relfval(end),output.relstep(end), ...
                    output.delta(end),output.rho(end));
        end
    end

end

% Display termination message.
if options.Display > 0
    ahref = '\n%s\n\n';
    x = round(linspace(0,output.iterations,min(500,output.iterations)));
    if length(bold) > 2
        ahref = sprintf(['\n<a href="matlab:semilogy(%s,%s);' ...
            'xlabel(''iteration'');legend(''fval'',' ...
            '''relfval'',''relstep'')">%%s</a>\n\n'],mat2str(x'), ...
            mat2str([output.fval(x+1)' [nan output.relfval(x(2:end))]' ...
                    [nan output.relstep(x(2:end))]'],3));
    end
    switch output.info
        case 1, fprintf(ahref,'Objective function tolerance reached.');
        case 2, fprintf(ahref,'Step size tolerance reached.');
        case 3, fprintf(ahref,'Maximum number of iterations reached.');
    end
end

end

function [z,offset] = deserialize(z,dim,offset)
    if iscell(dim)
        v = z;
        z = cell(size(dim));
        if nargin < 3, offset = 0; end
        for i = 1:numel(z)
            if iscell(dim{i})
                [z{i},offset] = deserialize(v,dim{i},offset);
            else
                n = prod(dim{i}(:));
                z{i} = reshape(v(offset+(1:n)),dim{i});
                offset = offset+n;
            end
        end
    elseif ~isempty(dim)
        z = reshape(z,dim);
    end
end

function z = serialize(z)
    if iscell(z)
        for i = find(cellfun(@iscell,z(:).'))
            z{i} = serialize(z{i});
        end
        s = cellfun(@numel,z(:)); o = [0; cumsum(s)];
        c = z; z = zeros(o(end),1);
        for i = 1:length(s), z(o(i)+(1:s(i))) = c{i}(:); end
    else
        z = z(:);
    end
end

function dim = structure(z)
    if iscell(z)
        dim = cellfun(@size,z,'UniformOutput',false);
        for i = find(cellfun(@iscell,z(:).'))
            dim{i} = structure(z{i});
        end
    else
        dim = size(z);
        if numel(z) == dim(1), dim = []; end
    end
end
