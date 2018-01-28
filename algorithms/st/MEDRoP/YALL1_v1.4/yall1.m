function [x, Out] = yall1(A, b, opts)
%
% A solver for L1-minimization models:
%
% min ||Wx||_{w,1}, st Ax = b
% min ||Wx||_{w,1} + (1/nu)||Ax - b||_1
% min ||Wx||_{w,1} + (1/2*rho)||Ax - b||_2^2
% min ||x||_{w,1}, st Ax = b                and x > = 0
% min ||x||_{w,1} + (1/nu)||Ax - b||_1,      st x > = 0
% min ||x||_{w,1} + (1/2*rho)||Ax - b||_2^2, st x > = 0
%
% where (A,b,x) can be complex or real 
% (but x must be real in the last 3 models)
%
% Copyright(c) 2009-2011 Yin Zhang, Junfeng Yang, Wotao Yin
%
% --- Input:
%     A --- either an m x n matrix or
%           a structure with 2 fields:
%           1) A.times: a function handle for A*x
%           2) A.trans: a function handle for A'*y
%     b --- an m-vector, real or complex
%  opts --- a structure with fields:
%           opts.tol   -- tolerance *** required ***
%           opts.nu    -- values > 0 for L1/L1 model
%           opts.rho   -- values > 0 for L1/L2 model
%           opts.basis -- sparsifying unitary basis W (W*W = I)
%                        a struct with 2 fields:
%                        1) times: a function handle for W*x
%                        2) trans: a function handle for W'*y
%           opts.nonneg  -- 1 for nonnegativity constraints
%           opts.nonorth -- 1 for A with non-orthonormal rows
%           see the User's Guide for all other options
%
% --- Output: 
%     x --- last iterate (hopefully an approximate solution)
%   Out --- a structure with fields:
%           Out.exit    --- exit information
%           Out.iter    --- #iterations taken
%           Out.cputime --- solver CPU time
%           Out.y       --- dual variable
%           Out.z       --- dual slack
%           .....       --- and some more
% --------------------------------------------------------------
% define linear operators
[A,At,b,opts] = linear_operators(A,b,opts);

m = length(b);
L1L1 = isfield(opts,'nu') && opts.nu > 0;
if L1L1 && isfield(opts,'weights')
    opts.weights = [opts.weights(:); ones(m,1)];
end

% parse options
posrho = isfield(opts,'rho')    && opts.rho > 0;
posdel = isfield(opts,'delta')  && opts.delta > 0;
posnu  = isfield(opts,'nu')     && opts.nu > 0;
nonneg = isfield(opts,'nonneg') && opts.nonneg == 1;
if isfield(opts,'x0'); x0 = opts.x0; else x0 = []; end 
if isfield(opts,'z0'); z0 = opts.z0; else z0 = []; end 

% check conflicts % modified by Junfeng
if posdel && posrho || posdel && posnu || posrho && posnu
    fprintf('Model parameter conflict! YALL1: set delta = 0 && nu = 0;\n');
    opts.delta = 0; posdel = false;
    opts.nu    = 0; posnu  = false;
end
prob = 'the basis pursuit problem';
if posrho, prob = 'the unconstrained L1L2 problem'; end
if posdel, prob = 'the constrained L1L2 problem';   end
if posnu,  prob = 'the unconstrained L1L1 problem'; end
%disp(['YALL1 is solving ', prob, '.']);

% check zero solution % modified by Junfeng
Atb = At(b);
bmax = norm(b,inf);
L2Unc_zsol = posrho && norm(Atb,inf) <= opts.rho;
L2Con_zsol = posdel && norm(b) <= opts.delta;
L1L1_zsol  = posnu  && bmax < opts.tol;
BP_zsol    = ~posrho && ~posdel && ~posnu && bmax < opts.tol;
zsol = L2Unc_zsol || L2Con_zsol || BP_zsol || L1L1_zsol;
if zsol  
    n = length(Atb);
    x = zeros(n,1); 
    Out.iter = 0;
    Out.cntAt = 1;
    Out.cntA = 0;
    Out.exit = 'Data b = 0';
    return; 
end
% ========================================================================

% scaling data and model parameters
b1 = b / bmax;
if posrho, opts.rho   = opts.rho / bmax; end
if posdel, opts.delta = opts.delta / bmax; end
if isfield(opts,'xs'), opts.xs = opts.xs/bmax; end
    
% solve the problem
t0 = cputime; 
[x1,Out] = yall1_solve(A, At, b1, x0, z0, opts);
Out.cputime = cputime - t0;

% restore solution x
x = x1 * bmax;
if L1L1; x = x(1:end-m); end
if isfield(opts,'basis')
    x = opts.basis.trans(x);
end
if nonneg; x = max(0,x); end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,At,b,opts] = linear_operators(A0, b0, opts)
%
% define linear operators A and At 
% (possibly modify RHS b if nu > 0)
%
b = b0;
if isnumeric(A0); 
    if size(A0,1) > size(A0,2); 
        error('A must have m < n');
    end
    A  = @(x) A0*x;
    At = @(y) (y'*A0)';
elseif isstruct(A0) && isfield(A0,'times') && isfield(A0,'trans');
    A  = A0.times;
    At = A0.trans;
elseif isa(A0,'function_handle')
    A  = @(x) A0(x,1);
    At = @(x) A0(x,2);
else
    error('A must be a matrix, a struct or a function handle');
end

% use sparsfying basis W
if isfield(opts,'basis')
    C = A; Ct = At; clear A At; 
    B  = opts.basis.times;
    Bt = opts.basis.trans;
    A  = @(x) C(Bt(x));
    At = @(y) B(Ct(y));
end

% solving L1-L1 model if nu > 0
if isfield(opts,'nu') && opts.nu > 0
    C = A; Ct = At; clear A At; 
    m = length(b0);
    nu = opts.nu; 
    t = 1/sqrt(1 + nu^2);
    A  = @(x) ( C(x(1:end-m)) + nu*x(end-m+1:end) )*t;
    At = @(y) [ Ct(y);  nu*y ]*t;
    b = b0*t;
end

if ~isfield(opts,'nonorth'); 
    opts.nonorth = check_orth(A,At,b); 
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nonorth = check_orth(A, At, b)
%
% check whether the rows of A are orthonormal
%
nonorth = 0;
s1 = randn(size(b));
s2 = A(At(s1));
err = norm(s1-s2)/norm(s1);
if err > 1.e-12; nonorth = 1; end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                   solver                %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, Out] = yall1_solve(A,At,b,x0,z0,opts)

% yall1_solve version 1.4 (July, 2011)
% Copyright(c) 2009-2011 Yin Zhang

%% initialization
m = length(b); bnrm = norm(b);
[tol,mu,maxit,print,nu,rho,delta, ... 
    w,nonneg,nonorth,gamma] = get_opts;
x = x0; z = z0;
if isempty(x0); x = At(b); end
n = length(x); 
if isempty(z0), z = zeros(n,1); end
if isfield(opts,'nonorth') && opts.nonorth > 0
    y = zeros(m,1); Aty = zeros(n,1); 
end
if print; fprintf('--- YALL1 v1.4 ---\n'); end
if print; iprint1(0); end;

mu_orig = mu;
rdmu = rho / mu;
rdmu1 = rdmu + 1;
bdmu = b / mu;
ddmu = delta / mu;

Out.cntA = 0; Out.cntAt = 0;
rel_gap = 0;  rel_rd  = 0;
rel_rp  = 0;  stop = 0;

%% main iterations
for iter = 1:maxit
    
    %% calculations
    xdmu = x / mu;
    if ~nonorth; % orthonormal A
        y = A(z - xdmu) + bdmu;
        if rho > 0;
            y = y / rdmu1;
        elseif delta > 0
            y = max(0, 1 - ddmu/norm(y))*y;
        end
        Aty = At(y);
    else     % non-orthonormal A
        ry = A(Aty - z + xdmu) - bdmu;
        if rho > 0; ry = ry + rdmu*y; end
        Atry = At(ry);
        denom = Atry'*Atry;
        if rho > 0, denom = denom + rdmu * (ry'*ry); end
        stp = real(ry'*ry)/(real(denom) + eps);
        Out.cntAt = Out.cntAt + 1;
        y = y - stp*ry;
        Aty = Aty - stp*Atry;
    end
    
    z = Aty + xdmu;
    z = proj2box(z,w,nonneg,nu,m);
    
    Out.cntA  = Out.cntA  + 1;
    Out.cntAt = Out.cntAt + 1;

    rd = Aty - z; xp = x;
    x = x + (gamma*mu) * rd;
        
    %% other chores
    if rem(iter,2) == 0, 
        check_stopping; update_mu; 
    end
    if print > 1; iprint2; end
    if stop; break; end 
    
end % main iterations

% output
Out.iter = iter;
Out.mu = [mu_orig mu];
Out.obj = [objp objd];
Out.y = y; Out.z = z;

if iter == maxit; Out.exit = 'Exit: maxiter'; end
if print; iprint1(1); end

%% nested functions
    function [tol,mu,maxit,print,nu,rho,delta, ...
             w,nonneg,nonorth,gamma] = get_opts
        % get or set options
        tol = opts.tol;
        mu = mean(abs(b)); 
        %mu = norm(b)/numel(b);
        maxit = 9999;
        print = 0;
        nu = 0;
        rho = eps;
        delta = 0;
        w = 1;
        nonneg = 0;
        nonorth = 0;
        gamma = 1.; % ADM parameter
        if isfield(opts,'mu');       mu = opts.mu;    end
        if isfield(opts,'maxit'); maxit = opts.maxit; end
        if isfield(opts,'print'); print = opts.print; end        
        if isfield(opts,'nu');       nu = opts.nu;    end
        if isfield(opts,'rho');     rho = opts.rho;   end
        if isfield(opts,'delta'); delta = opts.delta; end
        if isfield(opts,'weights'); w = opts.weights; end
        if isfield(opts,'nonneg');   nonneg = opts.nonneg;  end
        if isfield(opts,'nonorth'); nonorth = opts.nonorth; end
        if isfield(opts,'gamma');   gamma   = opts.gamma;   end
    end

    function z = proj2box(z,w,nonneg,nu,m)
        if nonneg
            z = min(w,real(z));
            if nu > 0 %L1L1 model
                z(end-m:end) = max(-1,z(end-m:end));
            end
        else
            z = z .* w ./ max(w,abs(z));
        end
    end

    function check_stopping
        q = 0.1; % q in [0,1)
        if delta > 0; q = 0; end
        % dual residual
        rdnrm = norm(rd); 
        rel_rd = rdnrm / norm(z);
        % duality gap
        objp = sum(abs(w.*x));
        objd = b'*y;
        if delta > 0, objd = objd - delta*norm(y); end
        if rho > 0
            rp = A(x) - b; 
            rpnrm = norm(rp); 
            Out.cntA = Out.cntA + 1;
            objp = objp + (0.5/rho)*rpnrm^2;
            objd = objd - (0.5*rho)*norm(y)^2;
        end
        rel_gap = abs(objd - objp)/abs(objp);
        
        % check relative change
        xrel_chg = norm(x-xp)/norm(x);
        if xrel_chg < tol*(1 - q)
            Out.exit = 'Exit: Stablized'; 
            stop = 1; return; 
        end
        
        % decide whether to go further
        if xrel_chg >= tol*(1 + q); return; end
        gap_small = rel_gap < tol;
        if ~gap_small; return; end
        d_feasible = rel_rd < tol;
        if ~d_feasible; return; end

        % check primal residual
        if rho == 0, 
            rp = A(x) - b; 
            rpnrm = norm(rp);
            Out.cntA = Out.cntA + 1; 
        end;    
        if rho > 0;
            p_feasible = true;
        elseif delta > 0
            p_feasible = rpnrm <= delta*(1 + tol);
        else
            p_feasible = rpnrm < tol*bnrm;
        end
        if p_feasible, stop = 1; Out.exit = 'Exit: Converged'; end        
    end

    function iprint1(mode)
        switch mode;
            case 0; % at the beginning
                rp = A(x) - b;
                rpnrm = norm(rp);
                fprintf(' norm( A*x0 - b ) = %6.2e\n',rpnrm);
            case 1; % at the end
                rp = A(x) - b;
                objp = sum(abs(w.*x));
                objd = b'*y;
                if rho > 0
                    objp = objp + (0.5/rho)*(rp'*rp);
                    objd = objd - (0.5*rho)*( y'*y );
                end
                if delta > 0; objd = objd - delta*norm(y); end
                dgap = abs(objd - objp);
                rel_gap = dgap / abs(objp);
                rdnrm = norm(rd);
                rel_rd = rdnrm / norm(z);
                rpnrm = norm(rp);
                rel_rp = rpnrm / bnrm;
                fprintf(' Rel_Gap   Rel_ResD  Rel_ResP\n');
                fprintf(' %8.2e  %8.2e  %8.2e\n',rel_gap,rel_rd,rel_rp);
        end
    end

    function iprint2
        rdnrm = norm(rd);
        rp = A(x) - b;
        rpnrm = norm(rp);
        objp = sum(abs(w.*x));
        objd = b'*y;
        if rho > 0
            objp = objp + (0.5/rho)*rpnrm^2;
            objd = objd - (0.5*rho)*(y'*y);
        end
        if delta > 0; 
            objd = objd - delta*norm(y); 
        end
        gap = abs(objd - objp);
        if ~rem(iter,50)
            fprintf('  Iter %4i:' ,iter);
            fprintf('  Rel_Gap  %6.2e',gap/objp);
            fprintf('  Rel_ResD %6.2e',rdnrm/norm(z));
            fprintf('  Rel_ResP %6.2e',norm(rp)/norm(b));
            fprintf('\n');
        end
        
        if ~isfield(opts,'xs'), return; end
        if isfield(opts,'nu') && opts.nu, return; end
        
        if iter <= 1, 
            Out.objd = []; Out.rd = []; 
            Out.objp = []; Out.rp = [];
            opts.xsnrm = norm(opts.xs);
            Out.relerr = [];
        end
        Out.objd = [Out.objd objd];
        Out.objp = [Out.objp objp];
        Out.rd = [Out.rd rdnrm];
        Out.rp = [Out.rp rpnrm];
        err = norm(x-opts.xs)/opts.xsnrm;
        Out.relerr = [Out.relerr err];
    end

    function update_mu    % added to v14
        mfrac = 0.1; big = 50; nup = 8;
        mu_min = mfrac^nup * mu_orig;
        do_update = rel_gap > big*rel_rd;
        do_update = do_update && mu > 1.1*mu_min;
        do_update = do_update && iter > 10;
        if ~do_update, return; end
        % do update
        mu = max(mfrac*mu,mu_min);
        rdmu = rho / mu; rdmu1 = rdmu + 1;
        bdmu = b / mu; ddmu = delta / mu;
        if print>1; fprintf('  -- mu updated\n'); end
    end

end
