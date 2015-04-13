function [z,output] = nlsb_gndl(F,dF,lb,ub,z0,options)
%NLSB_GNDL Bound-constrained NLS by projected Gauss-Newton dogleg TR.
%   [z,output] = nlsb_gndl(F,dF,lb,ub,z0) starts at z0 and attempts to find
%   a local minimizer of the real-valued function f(z), which is the
%   nonlinear least squares objective function f(z) := 0.5*(F(z)'*F(z))
%   subject to the constraints real(lb) <= real(z) <= real(ub) and
%   imag(lb) <= imag(z) <= imag(ub). The input variable z (and lb, ub) may
%   be a scalar, vector, matrix, tensor or even a (nested) cell array of
%   tensors and its contents may be real or complex. This method may be
%   applied in the following ways:
%
%   1. F is function of both z and conj(z).
%
%      Method 1: general medium-scale problems.
%      nlsb_gndl(F,dF,lb,ub,z0) where F(z) returns a column vector of
%      complex. Set dF equal to the string 'Jacobian-C' for automatic
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
%      nlsb_gndl(F,dF,lb,ub,z0) where F(z) returns a column vector of
%      complex residuals and dF is a structure containing:
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
%      nlsb_gndl(F,dF,lb,ub,z0) where F(z) returns a column vector of
%      complex residuals. Set dF equal to the string 'Jacobian' for
%      automatic numerical approximation of the Jacobian, respectively. Or,
%      supply the Jacobian manually with a structure dF containing:
%
%         dF.dz      - The function dF.dz(zk) should return the Jacobian
%                      dF(zk)/d(z^T), which is defined as the matrix in
%                      which the m-th row is equal to (dFm(zk)/dz)^T, where
%                      Fm is the m-th component of F.
%
%      Method 2: analytic large-scale problems.
%      nlsb_gndl(F,dF,lb,ub,z0) where F(z) returns a column vector of
%      complex residuals and dF is a structure containing:
%
%         dF.dzx     - The function dF.dzx(zk,x,'notransp') should return
%                      the matrix-vector product [dF(zk)/d(z^T)]*x and
%                      dF.dzx(zk,x,'transp') should return the matrix-
%                      vector product [dF(zk)/d(z^T)]'*x.
%
%      Method 3: analytic problems in a modest number of variables z and
%                large number of residuals F(z).
%      nlsb_gndl(f,dF,lb,ub,z0) where f(z) := 0.5*(F(z)'*F(z)) and dF is a
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
%      nlsb_gndl(f,dF,lb,ub,z0) where f(z) := 0.5*(F(z)'*F(z)) and dF is a
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
%      output.alpha        - The plane search step lengths in every
%                            iteration, if a plane search is selected.
%      output.cgiterations - The number of CG/LSQR iterations to compute
%                            the Gauss-Newton step in every iteration.
%                            (large-scale methods only).
%      output.cgrelres     - The relative residual norm of the computed
%                            Gauss-Newton step (large-scale methods only).
%      output.delta        - The trust region radius at every step attempt.
%      output.fval         - The value of the objective function f in every
%                            iteration.
%      output.info         - The circumstances under which the procedure
%                            terminated:
%                               1: Objective function tolerance reached.
%                               2: Step size tolerance reached.
%                               3: Maximum number of iterations reached.
%      output.infops       - The circumstances under which the plane search
%                            terminated in every iteration.
%      output.iterations   - The number of iterations.
%      output.relfval      - The difference in objective function value
%                            between every two successive iterates,
%                            relativeto its initial value.
%      output.relstep      - The step size relative to the norm of the 
%                            current iterate in every iteration.
%      output.rho          - The trustworthiness at every step attempt.
%
%   nlsb_gndl(F,dF,lb,ub,z0,options) may be used to set the following
%   options:
%
%      options.CGMaxIter = 15     - The maximum number of CG/LSQR
%                                   iterations for computing the
%                                   Gauss-Newton step (large-scale methods
%                                   only).
%      options.CGTol = 1e-6       - The tolerance for the CG/LSQR method to
%                                   compute the Gauss-Newton step
%                                   (large-scale methods only).
%      options.Delta =            - The initial trust region radius. If
%      0.3*max(1,norm(z0))          equal to NaN, the initial radius will
%                                   be equal to length of the first
%                                   Gauss-Newton step.
%      options.Display = 1        - Displays the objective function value,
%                                   its difference with the previous
%                                   iterate relative to the first iterate
%                                   and the relative step size each
%                                   options.Display iterations. Set to 0 to
%                                   disable.
%      options.JHasFullRank       - If set to true, the Gauss-Newton step
%      = false                      is computed as a least squares
%                                   solution, if possible. Otherwise, it is
%                                   computed using a more expensive
%                                   pseudo-inverse.
%      options.MaxIter = 200      - The maximum number of iterations.
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
%      options.TolFun = 1e-12     - The tolerance for output.relfval. Note
%                                   that because the objective function is
%                                   a squared norm, TolFun can be as small
%                                   as eps^2.
%      options.TolX = 1e-6        - The tolerance for output.relstep.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Unconstrained
%       optimization of real functions in complex variables", SIAM J. Opt.,
%       Vol. 22, No. 3, 2012, pp. 879-898.
%   [2] C.T. Kelley, "Iterative methods for optimization," SIAM Frontiers
%       in Applied Mathematics, No. 18, 1999.

% Check the objective function f, derivative dF and first iterate z0.
if ~isa(F,'function_handle')
    error('nlsb_gndl:F','The first argument must be a function.');
end
if ischar(dF)
    type = dF;
    if strcmp(type,'Jacobian-C'), fld = 'dzc'; else fld = 'dz'; end
    dF = struct(fld,@derivjac);
end
if ~isstruct(dF)
    error('nlsb_gndl:dF','Second argument not valid.');
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
        error('nlsb_gndl:dF', ...
             ['The structure dF should supply [dF.dzc] or ' ...
              '[dF.dzx and dF.dconjzx] or [dF.dz] or [dF.dzx] or ' ...
              '[dF.JHJ and dF.JHF] or [dF.JHJx and dF.JHF].']);
    end
end

% Define projection and reduction operators.
function z = proj(z)
    z = median([real(lb) real(z) real(ub)],2)+ ...
        median([imag(lb) imag(z) imag(ub)],2)*1i;
end
function H = red(H,A)
    H = bsxfun(@times,H,~A(1:size(H,1)));
    H = bsxfun(@times,H,~A(1:size(H,1)).');
    H(1:size(H,1)+1:end) = H(1:size(H,1)+1:end)+A(1:size(H,1));
end

% Evaluate the function value at z0.
lb = serialize(lb);
ub = serialize(ub);
dim = structure(z0);
z0 = proj(serialize(z0));
A = [real(z0)<=real(lb) | real(z0)>=real(ub) ...
     imag(z0)<=imag(lb) | imag(z0)>=imag(ub)];
z = deserialize(z0,dim);
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

% In the case 'F+dFdzx+dFdconjzx', compute a reduced version of J'*(J*x).
function y = JH_Jx(x)
    x1 = ~A(:,1).*x(1:end/2)+~A(:,2).*x(end/2+1:end)*1i;
    dFdzx = dF.dzx(z,x1,'notransp');
    dFdconjzconjx = dF.dconjzx(z,conj(x1),'notransp');
    y = real(dFdzx)+real(dFdconjzconjx)+ ...
        (imag(dFdzx)+imag(dFdconjzconjx))*1i;
	dFdzx = dF.dzx(z,y,'transp');
    dFdconjzconjx = dF.dconjzx(z,y,'transp');
    y = [~A(:,1).*(real(dFdzx)+real(dFdconjzconjx)); ...
         ~A(:,2).*(imag(dFdzx)-imag(dFdconjzconjx))];
    y(A(:)) = x(A(:));
end

% In the case 'F+dFdzx', compute a reduced version of dFdz'*(dFdz*x).
function y = dFdzH_dFdzx(x)
    x1 = ~A(:,1).*x(1:end/2)+~A(:,2).*x(end/2+1:end)*1i;
    y = dF.dzx(z,x1,'notransp');
    y = dF.dzx(z,y,'transp');
    y = [~A(:,1).*real(y);~A(:,2).*imag(y)];
    y(A(:)) = x(A(:));
end

% In the case 'f+JHJx+JHF', compute a reduced version of JHJ*x.
function y = JHJx(x)
    x1 = ~A(:,1).*x(1:end/2)+~A(:,2).*x(end/2+1:end)*1i;
    y = dF.JHJx(z,x1);
    y = [~A(:,1).*real(y);~A(:,2).*imag(y)];
    y(A(:)) = x(A(:));
end

% Modify the preconditioner, if available.
if isfield(dF,'M') && ~isempty(dF.M), dF.PC = @PC; else dF.PC = []; end
function x = PC(b)
    x = dF.M(z,b(1:end/2)+b(end/2+1:end)*1i);
    x = [real(x);imag(x)];
end

% Check the options structure.
if nargin < 6, options = struct; end
if ~isfield(options,'CGMaxIter'), options.CGMaxIter = 15; end
if ~isfield(options,'CGTol'), options.CGTol = 1e-6; end
if ~isfield(options,'Delta'), options.Delta = 0.3*max(1,norm(z0)); end
if ~isfield(options,'Display'), options.Display = 1; end
if ~isfield(options,'JHasFullRank'), options.JHasFullRank = false; end
if ~isfield(options,'MaxIter'), options.MaxIter = 200; end
if ~isfield(options,'PlaneSearch'), options.PlaneSearch = false; end
if ~isfield(options,'PlaneSearchOptions')
    options.PlaneSearchOptions = struct;
end
if ~isfield(options,'TolFun'), options.TolFun = 1e-12; end
if ~isfield(options,'TolX'), options.TolX = 1e-6; end

% Gauss-Newton with dogleg trust region.
output.alpha = [];
output.cgiterations = [];
output.cgrelres = [];
output.delta = options.Delta;
output.fval = fval;
output.info = false;
output.infops = [];
output.iterations = 0;
output.relfval = [];
output.relstep = [];
output.rho = [];
while ~output.info
    
    % Compute the (in)exact Gauss-Newton step pgn.
    switch method
        case 'F+dFdzc'
            % Compute the Gauss-Newton step pgn.
            dFdzc = dF.dzc(z);
            dFdz = dFdzc(:,1:end/2);
            dFdconjz = dFdzc(:,end/2+1:end);
            J = [real(dFdz)+real(dFdconjz),imag(dFdconjz)-imag(dFdz); ...
                 imag(dFdz)+imag(dFdconjz),real(dFdz)-real(dFdconjz)];
            grad = dFdz'*Fval+dFdconjz.'*conj(Fval);
            if options.JHasFullRank || issparse(J)
                pgn = red(J'*J,A)\[-real(grad);-imag(grad)];
            else
                pgn = pinv(red(J'*J,A))*[-real(grad);-imag(grad)];
            end
            pgn = pgn(1:end/2)+pgn(end/2+1:end)*1i;
            % Compute the Cauchy point pcp = -alpha*grad.
            gg = grad'*grad;
            gBg = dFdz*grad+dFdconjz*conj(grad);
            gBg = gBg'*gBg;
            alpha = gg/gBg;
        case 'F+dFdzx+dFdconjzx'
            % Compute the Cauchy point pcp = -alpha*grad.
            grad = dF.dzx(z,Fval,'transp')+ ...
                   conj(dF.dconjzx(z,Fval,'transp'));
            gg = grad'*grad;
            gBg = dF.dzx(z,grad,'notransp')+ ...
                  dF.dconjzx(z,conj(grad),'notransp');
            gBg = gBg'*gBg;
            alpha = gg/gBg;
            if ~isfinite(alpha), alpha = 1; end;
            % Compute the Gauss-Newton step pgn.
            [pgn,~,output.cgrelres(end+1),output.cgiterations(end+1)] = ...
                mpcg(@JH_Jx,[-real(grad);-imag(grad)], ...
                     options.CGTol,options.CGMaxIter,dF.PC,[], ...
                     -alpha*[real(grad);imag(grad)]);
            pgn = pgn(1:end/2)+pgn(end/2+1:end)*1i;
        case 'F+dFdz'
            % Compute the Gauss-Newton step pgn.
            dFdz = dF.dz(z);
            grad = dFdz'*Fval;
            JHJr = dFdz'*dFdz;
            gradr = -grad;
            JHJIsReal = isreal(JHJr);
            if ~JHJIsReal
                JHJr = [real(JHJr) -imag(JHJr); imag(JHJr) real(JHJr)];
                gradr = [-real(grad);-imag(grad)];
            end
            if options.JHasFullRank || issparse(dFdz)
                pgn = red(JHJr,A)\gradr;
            else
                pgn = pinv(red(JHJr,A))*gradr;
            end
            if ~JHJIsReal, pgn = pgn(1:end/2)+pgn(end/2+1:end)*1i; end
            % Compute the Cauchy point pcp = -alpha*grad.
            gg = grad'*grad;
            gBg = dFdz*grad;
            gBg = gBg'*gBg;
            alpha = gg/gBg;
        case 'F+dFdzx'
            % Compute the Cauchy point pcp = -alpha*grad.
            grad = dF.dzx(z,Fval,'transp');
            gg = grad'*grad;
            gBg = dF.dzx(z,grad,'notransp');
            gBg = gBg'*gBg;
            alpha = gg/gBg;
            if ~isfinite(alpha), alpha = 1; end;
            % Compute the Gauss-Newton step pgn.
            [pgn,~,output.cgrelres(end+1),output.cgiterations(end+1)] = ...
                mpcg(@dFdzH_dFdzx,-[real(grad);imag(grad)], ...
                     options.CGTol,options.CGMaxIter,dF.PC,[], ...
                     -alpha*[real(grad);imag(grad)]);
            pgn = pgn(1:end/2)+pgn(end/2+1:end)*1i;
        case 'f+JHJ+JHF'
            % Compute the Gauss-Newton step pgn.
            grad = serialize(dF.JHF(z));
            JHJ = dF.JHJ(z);
            JHJr = JHJ;
            gradr = -grad;
            JHJIsReal = isreal(JHJr);
            if ~JHJIsReal
                JHJr = [real(JHJr) -imag(JHJr); imag(JHJr) real(JHJr)];
                gradr = [-real(grad);-imag(grad)];
            end
            if options.JHasFullRank || issparse(JHJ)
                pgn = red(JHJr,A)\gradr;
            else
                pgn = pinv(red(JHJr,A))*gradr;
            end
            if ~JHJIsReal, pgn = pgn(1:end/2)+pgn(end/2+1:end)*1i; end
            % Compute the Cauchy point pcp = -alpha*grad.
            gg = grad'*grad;
            gBg = real(grad'*JHJ*grad);
            alpha = gg/gBg;
        case 'f+JHJx+JHF'
            % Compute the Cauchy point pcp = -alpha*grad.
            grad = serialize(dF.JHF(z));
            gg = grad'*grad;
            gBg = real(grad'*dF.JHJx(z,grad));
            alpha = gg/gBg;
            if ~isfinite(alpha), alpha = 1; end;
            % Compute the Gauss-Newton step pgn.
            [pgn,~,output.cgrelres(end+1),output.cgiterations(end+1)] = ...
                mpcg(@JHJx,-[real(grad);imag(grad)], ...
                     options.CGTol,options.CGMaxIter,dF.PC,[], ...
                     -alpha*[real(grad);imag(grad)]);
            pgn = pgn(1:end/2)+pgn(end/2+1:end)*1i;
    end
    
    % Project pgn and pcp.
    if ~all(isfinite(pgn)), pgn = -alpha*grad; end
    pgn = proj(z0+pgn)-z0;
    pcp = proj(z0-alpha*grad)-z0;
    gg = pcp'*pcp;
    
    % Plane search in the plane spanned by {pgn,pcp}.
    if ~all(isfinite(pgn)), pgn = -alpha*grad; end
    if isa(options.PlaneSearch,'function_handle')
        state = output; state.grad = grad;
        [alpha,outputps] = options.PlaneSearch( ...
            F,dF,z,deserialize(pgn,dim),deserialize(pcp,dim), ...
            state,options.PlaneSearchOptions);
        output.alpha(:,end+1) = alpha;
        if length(alpha) < 3, alpha(3) = 1; end
        p = alpha(1)*pgn+alpha(2)*pcp;
        z1 = deserialize(alpha(3)*(z0+p),dim);
        relstep = norm(p)/norm(z0); if isnan(relstep), relstep = 0; end
        if isfield(outputps,'fval'), fval = outputps.fval;
        else
            switch method
                case {'F+dFdzc','F+dFdzx+dFdconjzx','F+dFdz','F+dFdzx'}
                    Fval = F(z); fval = 0.5*sum(Fval(:)'*Fval(:));
                case {'f+JHJ+JHF','f+JHJx+JHF'}
                    fval = f(z);
            end
        end
        if isfield(outputps,'info')
            output.infops(end+1) = outputps.info;
        end
        rho = 1;
    else
        rho = -inf;
    end
    
    % Dogleg trust region.
    normpgn = norm(pgn);
    if isnan(output.delta(end)), output.delta(end) = max(1,normpgn); end
    while rho <= 0

        % Compute the dogleg step p.
        % Assume the projection did not alter the pgn or pcp too much,
        % computing dfval's would otherwise be quite expensive.
        delta = output.delta(end);
        if normpgn <= delta
            p = pgn;
        elseif sqrt(gg) >= delta
            p = delta/sqrt(gg)*pcp;
        else
            bma = pgn-pcp; bmabma = bma'*bma;
            c = real(pcp'*bma);
            if c <= 0
                beta = (-c+sqrt(c^2+bmabma*(delta^2-gg)))/bmabma;
            else
                beta = (delta^2-gg)/(c+sqrt(c^2+bmabma*(delta^2-gg)));
            end
            p = pcp+beta*bma;
        end
        
        % Estimate objective function improvement.
        % Because of projection, we can not avoid a matrix-vector product.
        switch method
            case 'F+dFdzc'
                dfval = dFdz*p+dFdconjz*conj(p);
                dfval = -real(p'*grad)-0.5*(dfval'*dfval);
            case 'F+dFdzx+dFdconjzx'
                dfval = dF.dzx(z,p,'notransp')+ ...
                        dF.dconjzx(z,conj(p),'notransp');
                dfval = -real(p'*grad)-0.5*(dfval'*dfval);
            case 'F+dFdz'
                dfval = dFdz*p;
                dfval = -real(p'*grad)-0.5*(dfval'*dfval);
            case 'F+dFdzx'
                dfval = dF.dzx(z,p,'notransp');
                dfval = -real(p'*grad)-0.5*(dfval'*dfval);
            case 'f+JHJ+JHF'
                dfval = -real(p'*grad)-0.5*real(p'*JHJ*p);
            case 'f+JHJx+JHF'
                dfval = -real(p'*grad)-0.5*real(p'*dF.JHJx(z,p));
        end

        % Compute the trustworthiness rho.
        if dfval > 0
            z1 = deserialize(z0+p,dim);
            switch method
                case {'F+dFdzc','F+dFdzx+dFdconjzx','F+dFdz','F+dFdzx'}
                    Fval = F(z1); Fval = Fval(:);
                    fval = 0.5*sum(Fval'*Fval);
                case {'f+JHJ+JHF','f+JHJx+JHF'}
                    fval = f(z1);
            end
            rho = (output.fval(end)-fval)/dfval;
            if isnan(rho), rho = -inf; end
            output.rho(end+1) = rho;
        end

        % Update trust region radius delta.
        if rho > 0.5
            output.delta(end+1) = max(delta,2*norm(p));
        else
            sigma = (1-0.25)/(1+exp(-14*(rho-0.25)))+0.25;
            if normpgn < sigma*delta && rho < 0
                e = ceil(log2(normpgn/delta)/log2(sigma));
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
            break;
        end

    end

    % Save current state.
    if rho > 0
        z = z1;
        z0 = serialize(z);
        A = [real(z0)<=real(lb) | real(z0)>=real(ub) ...
             imag(z0)<=imag(lb) | imag(z0)>=imag(ub)];
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
