function [z,output] = nls_lm(F,dF,z0,options)
%NLS_LM Nonlinear least squares by Levenberg-Marquardt.
%   [z,output] = nls_lm(F,dF,z0) starts at z0 and attempts to find a local
%   minimizer of the real-valued function f(z), which is the nonlinear
%   least squares objective function f(z) := 0.5*(F(z)'*F(z)). The input
%   variable z may be a scalar, vector, matrix, tensor or even a (nested)
%   cell array of tensors and its contents may be real or complex. This
%   method may be applied in the following ways:
%
%   1. F is function of both z and conj(z).
%
%      Method 1: general medium-scale problems.
%      nls_lm(F,dF,z0) where F(z) returns a column vector of complex
%      residuals. Set dF equal to the string 'Jacobian-C' for automatic
%      numerical approximation of the complex Jacobian, or supply the
%      complex Jacobian manually with a structure dF containing:
%
%         dF.dzc     - The function dF.dzc(zk) should return the complex
%                      Jacobian [dF(zk)/d(z^T) dF(zk)/d(conj(z)^T)], which
%                      is defined as the matrix in which the m-th row is
%                      equal to [(dFm(zk)/dz); (dFm(zk)/d(conj(z)))]^T,
%                      where Fm is the m-th component of F.
%
%      Method 2: general large-scale problems.
%      nls_lm(F,dF,z0) where F(z) returns a column vector of complex
%      residuals and dF is a structure containing:
%
%         dF.dzx     - The function dF.dzx(zk,x,'notransp') should return
%                      the matrix-vector product [dF(zk)/d(z^T)]*x and
%                      dF.dzx(zk,x,'transp') should return the
%                      matrix-vector product [dF(zk)/d(z^T)]'*x.
%         dF.dconjzx - The function dF.dconjzx(zk,x,'notransp') should
%                      return the matrix-vector product
%                      [dF(zk)/d(conj(z)^T)]*x and
%                      dF.dconjzx(zk,x,'transp') should return the matrix-
%                      vector product [dF(zk)/d(conj(z)^T)]'*x.
%
%   2. F is function only of z.
%
%      Method 1: analytic medium-scale problems.
%      nls_lm(F,dF,z0) where F(z) returns a column vector of complex
%      residuals. Set dF equal to the string 'Jacobian' for automatic
%      numerical approximation of the Jacobian, respectively. Or, supply
%      the Jacobian manually with a structure dF containing:
%
%         dF.dz      - The function dF.dz(zk) should return the Jacobian
%                      dF(zk)/d(z^T), which is defined as the matrix in
%                      which the m-th row is equal to (dFm(zk)/dz)^T, where
%                      Fm is the m-th component of F.
%
%      Method 2: analytic large-scale problems.
%      nls_lm(F,dF,z0) where F(z) returns a column vector of complex
%      residuals and dF is a structure containing:
%
%         dF.dzx     - The function dF.dzx(zk,x,'notransp') should return
%                      the matrix-vector product [dF(zk)/d(z^T)]*x and
%                      dF.dzx(zk,x,'transp') should return the matrix-
%                      vector product [dF(zk)/d(z^T)]'*x.
%
%      Method 3: analytic problems in a modest number of variables z and
%                large number of residuals F(z).
%      nls_lm(f,dF,z0) where f(z) := 0.5*(F(z)'*F(z)) and dF is a
%      structure containing:
%
%         dF.JHF     - The function dF.JHF(zk) should return 
%                      [dF(zk)/d(z^T)]'*F(zk), which is also equal to 
%                      2*df(zk)/d(conj(z)) = 2*conj(df(zk)/d(z)) if z is
%                      complex, or equal to df(xk)/dx if it is real.
%         dF.JHJ     - The function dF.JHF(zk) should return the Gramian
%                      [dF(zk)/d(z^T)]'*[dF(zk)/d(z^T)].
%      
%      Method 4: analytic problems in a large number of variables z and
%                large number of residuals F(z).
%      nls_lm(f,dF,z0) where f(z) := 0.5*(F(z)'*F(z)) and dF is a
%      structure containing:
%
%         dF.JHF     - The function dF.JHF(zk) should return 
%                      [dF(zk)/d(z^T)]'*F(zk), which is also equal to 
%                      2*df(zk)/d(conj(z)) = 2*conj(df(zk)/d(z)) if z is
%                      complex, or equal to df(xk)/dx if it is real.
%         dF.JHJx    - The function dF.JHF(zk,x) should return the matrix-
%                      vector product ([dF(zk)/d(z^T)]'*[dF(zk)/d(z^T)])*x.
%
%   The structure output returns additional information:
%
%      output.cgiterations - The number of CG/LSQR iterations to compute 
%                            the Levenberg-Marquardt step in every
%                            iteration (large-scale methods only).
%      output.cgrelres     - The relative residual norm of the computed 
%                            Levenberg-Marquardt step (large-scale methods
%                            only).
%      output.fval         - The value of the objective function f in every
%                            iteration.
%      output.info         - The circumstances under which the procedure
%                            terminated:
%                               1: Objective function tolerance reached.
%                               2: Step size tolerance reached.
%                               3: Maximum number of iterations reached.
%      output.iterations   - The number of iterations.
%      output.lambda       - The regularization factor lambda of the
%                            modified objective function 0.5*(F(z)'*F(z)+
%                            lambda*z'*z) at every step attempt.
%      output.relfval      - The difference in objective function value
%                            between every two successive iterates,
%                            relativeto its initial value.
%      output.relstep      - The step size relative to the norm of the 
%                            current iterate in every iteration.
%      output.rho          - The trustworthiness at every step attempt.
%
%   nls_lm(F,dF,z0,options) may be used to set the following options:
%
%      options.CGMaxIter = 15 - The maximum number of CG/LSQR iterations
%                               for computing the Levenberg-Marquardt step 
%                               (large-scale methods only).
%      options.CGTol = 1e-6   - The tolerance for the CG/LSQR method to
%                               compute the Levenberg-Marquardt step
%                               (large-scale methods only).
%      options.Display = 1    - Displays the objective function value, its
%                               difference with the previous iterate
%                               relative to the first iterate and the
%                               relative step size each options.Display
%                               iterations. Set to 0 to disable.
%      options.Lambda = 1     - The initial regularization factor lambda
%                               of the modified objective function
%                               0.5*(F(z)'*F(z)+lambda*z'*z).
%      options.MaxIter = 200  - The maximum number of iterations.
%      options.TolFun = 1e-12 - The tolerance for output.relfval. Note that
%                               because the objective function is a squared
%                               norm, TolFun can be as small as eps^2.
%      options.TolX = 1e-6    - The tolerance for output.relstep.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Unconstrained
%       optimization of real functions in complex variables", SIAM J. Opt.,
%       Vol. 22, No. 3, 2012, pp. 879-898.

% Check the objective function f, derivative dF and first iterate z0.
if ~isa(F,'function_handle')
    error('nls_lm:F','The first argument must be a function.');
end
if ischar(dF)
    type = dF;
    if strcmp(type,'Jacobian-C'), fld = 'dzc'; else fld = 'dz'; end
    dF = struct(fld,@derivjac);
end
if ~isstruct(dF)
    error('nls_lm:dF','Second argument not valid.');
else
    if isfield(dF,'dzc')
        method = 'F+dFdzc';
    elseif isfield(dF,'dzx') && isfield(dF,'dconjzx')
        method = 'F+dFdzx+dFdconjzx';
    elseif isfield(dF,'dz')
        method = 'F+dFdz';
    elseif isfield(dF,'dzx')
        method = 'F+dFdzx';
    elseif isfield(dF,'JHJ')  && isfield(dF,'JHF')
        method = 'f+JHJ+JHF';
        f = F;
    elseif isfield(dF,'JHJx') && isfield(dF,'JHF')
        method = 'f+JHJx+JHF';
        f = F;
    else
        error('nls_lm:dF', ...
             ['The structure dF should supply [dF.dzc] or ' ...
              '[dF.dzx and dF.dconjzx] or [dF.dz] or [dF.dzx] or ' ...
              '[dF.JHJ and dF.JHF] or [dF.JHJx and dF.JHF].']);
    end
end

% Evaluate the function value at z0.
dim = structure(z0);
z = z0;
z0 = serialize(z0);
n = length(z0);
switch method
    case {'F+dFdzc','F+dFdzx+dFdconjzx','F+dFdz','F+dFdzx'}
        Fval = F(z); Fval = Fval(:);
        fval = 0.5*sum(Fval'*Fval);
    case {'f+JHJ+JHF','f+JHJx+JHF'}
        fval = f(z);
end

% Numerical approximaton of complex derivatives.
function J = derivjac(zk)
    J = deriv(F,zk,Fval,type);
end

% In the case 'F+dFdzx+dFdconjzx', convert J(')*x to the real domain.
function y = Jx(x,transp)
    if strcmp(transp,'notransp')
        x = x(1:end/2)+x(end/2+1:end)*1i;
        dFdzx = dF.dzx(z,x,transp);
        dFdconjzconjx = dF.dconjzx(z,conj(x),transp);
        y = [real(dFdzx)+real(dFdconjzconjx); ...
             imag(dFdzx)+imag(dFdconjzconjx); ...
             sqrt(lambda)*real(x); ...
             sqrt(lambda)*imag(x)];
    else
        Jx = x(1:end/4)+x(end/4+1:end/2)*1i;
        lx = x(end/2+1:end*3/4)+x(end*3/4+1:end)*1i;
        dFdzx = dF.dzx(z,Jx,transp);
        dFdconjzx = dF.dconjzx(z,Jx,transp);
        y = [real(dFdzx)+real(dFdconjzx)+sqrt(lambda)*real(lx); ...
             imag(dFdzx)-imag(dFdconjzx)+sqrt(lambda)*imag(lx)];
    end
end

% In the case 'F+dFdzx', compute [dFdz;sqrt(lambda)*pad_eye](')*x.
function y = dFdzx(x,transp)
    if strcmp(transp,'notransp')
        y = [dF.dzx(z,x,transp);sqrt(lambda)*x];
    else
        Jx = x(1:end/2);
        lx = x(end/2+1:end);
        y = dF.dzx(z,Jx,transp)+sqrt(lambda)*lx;
    end
end

% In the case 'f+JHJx+JHF', compute (JHJ+lambda*pad_eye)*x.
function y = JHJx(x)
    y = dF.JHJx(z,x)+lambda*x;
end

% Modify the preconditioner, if available.
if isfield(dF,'M') && ~isempty(dF.M), dF.PC = @PC; else dF.PC = []; end
function x = PC(b)
    x = dF.M(z,b);
end

% Check the options structure.
if nargin < 4, options = struct; end
if ~isfield(options,'CGMaxIter'), options.CGMaxIter = 15; end
if ~isfield(options,'CGTol'), options.CGTol = 1e-6; end
if ~isfield(options,'Display'), options.Display = 1; end
if ~isfield(options,'Lambda'), options.Lambda = 1; end
if ~isfield(options,'MaxIter'), options.MaxIter = 200; end
if ~isfield(options,'TolFun'), options.TolFun = 1e-12; end
if ~isfield(options,'TolX'), options.TolX = 1e-6; end

% Levenberg-Marquardt.
output.cgiterations = [];
output.cgrelres = [];
output.fval = fval;
output.info = false;
output.iterations = 0;
output.lambda = options.Lambda;
output.relfval = [];
output.relstep = [];
output.rho = [];
while ~output.info
    
    % Evaluate first-order information.
    switch method
        case 'F+dFdzc'
            dFdzc = dF.dzc(z);
            dFdz = dFdzc(:,1:end/2);
            dFdconjz = dFdzc(:,end/2+1:end);
            if output.iterations == 0
                pad_vec = zeros(2*n,1);
                pad_eye = eye(2*n);
                if issparse(dFdz) || n < size(dFdz,1)
                    pad_eye = speye(2*n);
                end
            end
            J = [real(dFdz)+real(dFdconjz),imag(dFdconjz)-imag(dFdz); ...
                 imag(dFdz)+imag(dFdconjz),real(dFdz)-real(dFdconjz)];
            grad = dFdz'*Fval+dFdconjz.'*conj(Fval);
            b = [-real(Fval);-imag(Fval);pad_vec];
        case 'F+dFdzx+dFdconjzx'
            if output.iterations == 0, pad_vec = zeros(2*n,1); end
            grad = dF.dzx(z,Fval,'transp')+ ...
                   conj(dF.dconjzx(z,Fval,'transp'));
            b = [-real(Fval);-imag(Fval);pad_vec];
        case 'F+dFdz'
            dFdz = dF.dz(z);
            if output.iterations == 0
                pad_vec = zeros(n,1);
                pad_eye = eye(n);
                if issparse(dFdz) || n < size(dFdz,1)
                    pad_eye = speye(n);
                end
            end
            grad = dFdz'*Fval;
            b = [-Fval;pad_vec];
        case 'F+dFdzx'
            if output.iterations == 0, pad_vec = zeros(n,1); end
            grad = dF.dzx(z,Fval,'transp');
            b = [-Fval;pad_vec];
        case 'f+JHJ+JHF'
            JHJ = dF.JHJ(z);
            if output.iterations == 0
                pad_eye = eye(n);
                if issparse(JHJ), pad_eye = speye(n); end
            end
            grad = serialize(dF.JHF(z));
        case 'f+JHJx+JHF'
            grad = serialize(dF.JHF(z));
    end

    % Levenberg-Marquardt.
    rho = -inf;
    while rho <= 0
  
        % Compute the (in)exact Levenberg-Marquardt step plm.
        lambda = output.lambda(end);
        switch method
            case 'F+dFdzc'
                if n >= size(dFdz,1)
                    plm = [J;sqrt(lambda)*pad_eye]\b;
                else
                    plm = (J'*J+lambda*pad_eye)\[-real(grad);-imag(grad)];
                end
                plm = plm(1:end/2)+plm(end/2+1:end)*1i;
            case 'F+dFdzx+dFdconjzx'
                [plm,~,output.cgrelres(end+1), ...
                       output.cgiterations(end+1)] = ...
                    lsqr(@Jx,b,options.CGTol,options.CGMaxIter,dF.PC);
                plm = plm(1:end/2)+plm(end/2+1:end)*1i;
            case 'F+dFdz'
                if n >= size(dFdz,1)
                    plm = [dFdz;sqrt(lambda)*pad_eye]\b;
                else
                    plm = (dFdz'*dFdz+lambda*pad_eye)\(-grad);
                end
            case 'F+dFdzx'
                [plm,~,output.cgrelres(end+1), ...
                       output.cgiterations(end+1)] = ...
                    lsqr(@dFdzx,b,options.CGTol,options.CGMaxIter,dF.PC);
            case 'f+JHJ+JHF'
                plm = (JHJ+lambda*pad_eye)\(-grad);
            case 'f+JHJx+JHF'
                [plm,~,output.cgrelres(end+1), ...
                       output.cgiterations(end+1)] = ...
                   mpcg(@JHJx,-grad,options.CGTol,options.CGMaxIter,dF.PC);
        end
        
        % Compute the trustworthiness rho.
        if any(isnan(plm)), plm = zeros(size(plm)); end
        dfval = 0.5*abs(lambda*(plm'*plm)-real(grad'*plm));
        if dfval > 0
            z = deserialize(z0+plm,dim);
            switch method
                case {'F+dFdzc','F+dFdzx+dFdconjzx','F+dFdz','F+dFdzx'}
                    Fval = F(z); Fval = Fval(:);
                    fval = 0.5*sum(Fval'*Fval);
                case {'f+JHJ+JHF','f+JHJx+JHF'}
                    fval = f(z);
            end
            rho = (output.fval(end)-fval)/dfval;
            if isnan(rho), rho = -inf; end
            output.rho(end+1) = rho;
        end
        
        % Update lambda.
        if rho > 0
            output.lambda(end+1) = lambda*max(1/3,1-(2*rho-1)^3);
            nu = 2;
        else
            if output.iterations == 0, nu = 2; end
            output.lambda(end+1) = lambda*nu;
            nu = 2*nu;
        end
        
        % Check for convergence.
        relstep = norm(plm)/norm(z0); if isnan(relstep), relstep = 0; end
        if rho <= 0 && relstep <= options.TolX
            output.rho(end+1) = rho;
            fval = output.fval(end);
            z = deserialize(z0,dim);
            break;
        end
    
    end
    
    % Save current state.
    if rho > 0
        z0 = z0+plm;
    end
    
    % Update the output structure.
    output.fval(end+1) = fval;
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
            fprintf('%10s%s','',sprintf(bold,'lambda'));
            fprintf('%7s%s','',sprintf(bold,'rho'));
            fprintf('\n%21s%9s = %4.e %6s = %4.e\n\n','=1/2*norm(F)^2', ...
                    'TolFun',options.TolFun,'TolX',options.TolX);
        end
        if output.iterations == 1
            fprintf('%4i: % 14.8e |\n',0,output.fval(1));
        end
        fprintf('%4i: % 14.8e | %14.8e | %14.8e | %10.4e | %10.4e\n', ...
                output.iterations,output.fval(end), ...
                output.relfval(end),output.relstep(end), ...
                output.lambda(end),output.rho(end));
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
