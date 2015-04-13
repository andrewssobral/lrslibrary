function [z,output] = minf_ncg(f,g,z0,options)
%MINF_NCG Minimize a function by nonlinear conjugate gradient.
%   [z,output] = minf_ncg(f,g,z0) starts at z0 and attempts to find a local
%   minimizer of the real-valued function f(z). The input variables z may
%   be a scalar, vector, matrix, tensor or even a (nested) cell array of
%   tensors and its contents may be real or complex.
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
%      output.alpha      - The line search step length in every iteration.
%      output.fevals     - The total number of function/gradient calls.
%      output.fval       - The value of the objective function f in every
%                          iteration.
%      output.info       - The circumstances under which the procedure
%                          terminated:
%                             1: Objective function tolerance reached.
%                             2: Step size tolerance reached.
%                             3: Maximum number of iterations reached.
%      output.infols     - The circumstances under which the line search
%                          terminated in every iteration.
%      output.iterations - The number of iterations.
%      output.relfval    - The difference in objective function value
%                          between every two successive iterates, relative
%                          to its initial value.
%      output.relstep    - The step size relative to the norm of the 
%                          current iterate in every iteration.
%
%   minf_ncg(f,g,z0,options) may be used to set the following options:
%
%      options.Beta =            - The conjugate gradient update parameter.
%      ['HS',{'HSm'},'PR',...      Options are 'HS' (Hestenes-Stiefel),
%       'PRm','FR,'SD']            'HSm' (modified Hestenes-Stiefel), 'PR'
%                                  (Polak-Ribiere), 'PRm' (modified Polak-
%                                  Ribiere), 'FR' (Fletcher-Reeves) and
%                                  'SD' (steepest descent).
%      options.Display = 10      - Displays the objective function value,
%                                  its difference with the previous iterate
%                                  relative to the first iterate and the
%                                  relative step size each options.Display
%                                  iterations. Set to 0 to disable.
%      options.LineSearch        - The line search used to minimize the
%      = @ls_mt                    objective function in the quasi-Newton
%                                  descent direction.
%      options.LineSearchOptions - The options structure passed to the line
%                                  search routine.
%      options.MaxIter = 500     - The maximum number of iterations.
%      options.TolFun = 1e-6     - The tolerance for output.relfval.
%      options.TolX = 1e-8       - The tolerance for output.relstep.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Unconstrained
%       optimization of real functions in complex variables", SIAM J. Opt.,
%       Vol. 22, No. 3, 2012, pp. 879-898.

% Evaluate the objective function and gradient.
dim = structure(z0);
fval = f(z0);
if ~isa(g,'function_handle') && isempty(g)
    grad = serialize(deriv(f,z0,fval));
else
    grad = serialize(g(z0));
end
z = z0;

% Check the options structure.
isfunc = @(f)isa(f,'function_handle');
if nargin < 4, options = struct; end
if ~isfield(options,'alpha')
    if ~isreal(grad), options.LineSearchOptions.alpha = 0.5;
    else options.alpha = 1; end
end
if ~isfield(options,'Beta'), options.Beta = 'HSm'; end
if ~isfield(options,'Display'), options.Display = 10; end
if ~isfield(options,'LineSearch') || ~isfunc(options.LineSearch)
    options.LineSearch = @ls_mt;
end
if ~isfield(options,'LineSearchOptions')
    options.LineSearchOptions = struct;
end
if ~isfield(options.LineSearchOptions,'alpha')
    options.LineSearchOptions.alpha = 1e-4;
end
if ~isfield(options.LineSearchOptions,'c2')
    options.LineSearchOptions.c2 = 0.1;
end
if ~isfield(options,'MaxIter'), options.MaxIter = 500; end
if ~isfield(options,'TolFun'), options.TolFun = 1e-6; end
if ~isfield(options,'TolX'), options.TolX = 1e-8; end

% Select the conjugate gradient update parameter.
switch options.Beta
    % Hestenes-Stiefel.
    case 'HS',  beta = @(g,g1,p)real((g-g1)'*g)/real((g-g1)'*p);
    case 'HSm', beta = @(g,g1,p)max(0,real((g-g1)'*g)/real((g-g1)'*p));
    % Polak-Ribiere (equivalent to HS for exact line search).
    case 'PR',  beta = @(g,g1,~)real((g-g1)'*g)/(g1'*g1);
    case 'PRm', beta = @(g,g1,~)max(0,real((g-g1)'*g)/(g1'*g1));
    % Fletcher-Reeves (equivalent to PR for quadratic objective functions).
    case 'FR',  beta = @(g,g1,~)g'*g/(g1'*g1);
    % Steepest descent.
    case 'SD',  beta = @(~,~,~)0;
    otherwise,  beta = options.Beta;
end

% Nonlinear conjugate gradient.
output.alpha = [];
output.fevals = 1;
output.fval = fval;
output.info = false;
output.infols = [];
output.iterations = 0;
output.relfval = [];
output.relstep = [];
while ~output.info

    % Compute the descent direction pncg.
    if output.iterations == 0
        pncg = -grad;
    else
        pncg = -grad+beta(grad,grad1,pncg)*pncg;
    end
    
    % Scale the initial step length.
    gp = real(grad'*pncg);
    if gp > 0
        pncg = -grad;
        options.LineSearchOptions.alpha = output.alpha(1);
    elseif output.iterations > 0 && g1p1 < 0
        options.LineSearchOptions.alpha = ...
            output.alpha(end)*max(min(2,gp/g1p1),0.05);
    end
    g1p1 = gp;
    
    % Minimize f along z+alpha*pncg.
    state = output; state.grad = grad;
    [alpha,outputls] = options.LineSearch( ...
        f,g,z,deserialize(pncg,dim),state,options.LineSearchOptions);
    output.alpha(:,end+1) = alpha;
    if length(alpha) < 2, alpha(2) = 1; end
    
    % Update iterate.
    z = deserialize(alpha(2)*(serialize(z)+alpha(1)*pncg),dim);
    
    % Update gradient and Hessian approximation.
    grad1 = grad;
    if isfield(outputls,'grad')
        grad = outputls.grad;
    else
        if ~isa(g,'function_handle') && isempty(g)
            grad = serialize(deriv(f,z,output.fval(end)));
        else
            grad = serialize(g(z));
        end
    end
    
    % Update the output structure.
    if isfield(outputls,'fevals')
        output.fevals = output.fevals+outputls.fevals;
    end
    if isfield(outputls,'fval')
        output.fval(end+1) = outputls.fval;
    else
        output.fval(end+1) = f(z);
    end
    if isfield(outputls,'info')
        output.infols(end+1) = outputls.info;
    end
    output.iterations = output.iterations+1;
    output.relfval(end+1) = ...
        abs(diff(output.fval(end:-1:end-1)))/abs(output.fval(1));
    output.relstep(end+1) = ...
        norm(output.alpha(end)*pncg)/norm(serialize(z));
    if isnan(output.relstep(end)), output.relstep(end) = 0; end
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
            fprintf('%10s%s','',sprintf(bold,'alpha'));
            fprintf('\n%30s = %4.e %6s = %4.e\n\n', ...
                    'TolFun',options.TolFun,'TolX',options.TolX);
        end
        if output.iterations == 1
            fprintf('%4i: % 14.8e |\n',0,output.fval(1));
        end
        stralpha = [repmat('%10.4e ',1,size(output.alpha,1)) '\n'];
        fprintf(['%4i: % 14.8e | %14.8e | %14.8e | ' stralpha], ...
                output.iterations,output.fval(end), ...
                output.relfval(end),output.relstep(end), ...
                abs(output.alpha(:,end)));
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
