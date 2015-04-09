function [U,S,V,Out] = LMSVD(A,r,opts)
%
% LMSVD  Limited Memory Block Krylov Subspace Optimization for
%        Computing Principal Singular Value Decompositions
% 
%     Input:
%        A --- either a numeric matrix or a struct with 4 fields:
%            A.size  -- [m, n]
%            A.times -- a function handle for A*x
%            A.trans -- a function handle for A'*x
%            A.param -- parameters used by A.times and A.trans
%                       (assign [] if none)
%            r --- number of leading singular triplets
%       opts --- option structure with fields (default in []):
%            tol     -- tolerance [1.e-8]
%            maxit   -- maximal number of iteration [300]
%            memo    -- number of block subspaces
%                       [default is set from relations between r and m,n]
%            gvk     -- number of addional guard vectors [10]
%            initY   -- initial guess of a n by r matrix [randn(n,r)]                      
%            idisp   -- detailed information display option [0]
%     Output:
%        U,S,V --- principal SVD of A with the r singular triplets
%        Out   --- output information
% 
%     Copyright 2012-2014. Versions 1.0 - 1.1.
%     Written by Xin Liu, Yin Zhang and Zaiwen Wen. July, 2010.
%     Revised: March, 2012.
%     Revised:  June, 2014.

if nargin < 2; r = 6;  end
[m,n] = check_matrix(A);
if r > min(m,n)/2;
    warning(TooManySVSrequested,'r > min(m,n)/2');
end

% set parameters
mainargin = nargin;
set_param;

% initialize
if isnumeric(A)
    tA1 = tic; X = A*Y;         tA1 = toc(tA1);
    tqr = tic; [X,R] = qr(X,0); tqr = toc(tqr);
    tA2 = tic; Y = (X'*A)';     tA2 = toc(tA2);
else
    tA1 = tic; X = feval(A.times,Y,A.param); tA1 = toc(tA1);
    tqr = tic; [X,R] = qr(X,0);              tqr = toc(tqr);
    tA2 = tic; Y = feval(A.trans,X,A.param); tA2 = toc(tA2);
end

% bound memo
tAs = max(4*eps,(tA1 + tA2)/2);
tqr = max(  eps,tqr);
memb = ceil(tAs/tqr) + 1;
memo = max(0,min(memo,memb));

% call solver
[X,Y,Out] = lm_lbo(A,X,Y,r,tol,maxit,memo,idisp);

% generate svd
[U,S,V] = get_svd(X,Y);
Out.svk = diag(S);
Out.hrvs(end,:) = Out.svk(r:-1:1).^2;
%
% output principal SVD
U = U(:,1:r);
V = V(:,1:r);
S = S(1:r,1:r);
%%%% end of the main program %%%%


%% %%%%%%% nested functions %%%%%%% %%
    function set_param
        
        tol = 1e-8;
        maxit = 300;
        idisp = 0;
        mn = min(m,n);
        if r <= mn*0.02;
            memo = 5;
        elseif r <= mn*0.03;
            memo = 4;
        else
            memo = 3;
        end
        
        if isfield(opts,  'gvk'); tau = opts.gvk; else tau = 10; end
        % working size
        k = min([2*r,r+tau,m,n]);
        % initial guess
        Y = randn(n,k);
        
        if mainargin < 3; return; end  
        if isfield(opts,  'tol');     tol = opts.tol;   end
        if isfield(opts,'maxit');   maxit = opts.maxit; end
        if isfield(opts, 'memo');    memo = opts.memo;  end
        if isfield(opts,'idisp');   idisp = opts.idisp; end
        if isfield(opts,'initY');       Y = opts.initY(:,1:k); end
        
    end % set_param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [U,S,V] = get_svd(X,Y)
        method = 2;
        switch method
            case 1;
                [V,S,W] = svd(Y,0);
                U = X*W;
            case 2;
                [V,R] = qr(Y,0);
                [W,S,Z] = svd(R');
                U = X*W; V = V*Z;
        end
    end % get svd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [m,n] = check_matrix(A)
        if isnumeric(A)
            [m,n] = size(A);
        else % checking struct A
            if ~isstruct(A)
                error('A must be either numeric or struct');
            end
            if ~isfield(A,'times'); error('A.times missing'); end
            if ~isfield(A,'trans'); error('A.trans missing'); end
            if ~isfield(A,'size');  error('A.size  missing'); end
            if ~isa(A.times,'function_handle')
                error('A.times is not a function handle');
            end
            if ~isa(A.trans,'function_handle')
                error('A.trans is not a function handle');
            end
            m = A.size(1); n = A.size(2);
        end
    end % check A

%%%%%%%%%%%%
end % main %
%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% solver %%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,Out] = lm_lbo(A,X,Y,r,tol,maxit,memo,idisp)
%
% This code solves
%       min ||XY'-A||_F, s.t. X'^X = I,
% using a limited memory look-back optimization (LBO)
% acceleration.  The problem is equivalent to
%       max ||A'*X||_F,  s.t. X'*X = I.
%
% Input required:
%       A --- an (m by n) matrix or a struct
%       X --- an (m by k) matrix so that X'*X=I
%       Y --- an (n by k) matrix so that Y = A'*X
%

m = size(X,1); n = size(Y,1);
mn = min(m,n); k = size(Y,2);
if k < r; error('working size too small'); end

Xm = zeros(m,(1+memo)*k);
Ym = zeros(n,(1+memo)*k);
Xm(:,k+1:2*k) = X;
Ym(:,k+1:2*k) = Y;
Lm = k;

rvr = zeros(r,1);
chg_rvr = 1;
chgv = zeros(maxit,1);
xtrm = zeros(maxit,1);
hrvs = zeros(maxit,r);
kktc = zeros(maxit,1);
disp_str = 'iter %3i: memo used %i, chg_rvr %8.4e\n';

% set tolerance for terminating criterion
set_tolerance;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iter = 1:maxit
    SX = X;
    %% subspace iteration
    if isnumeric(A)
        AY = A*Y;
        [X,~] = qr(AY,0);
        Y = (X'*A)';
    else
        AY = feval(A.times,Y,A.param);
        [X,~] = qr(AY,0);
        Y = feval(A.trans,X,A.param);
    end % ---------------------------
    
    %% calculating terminating rules
    if Lm == 0 || iter <= 3
        SYTY = SX'*AY; SYTY = 0.5*(SYTY+SYTY');
        [tU,tE] = eig(SYTY);
        rvr0 = rvr;
        rv_sort = sort(diag(tE),'ascend');
        rvr = rv_sort(end-r+1:end);
        chg_rvr = norm(rvr0-rvr)/norm(rvr);
        hrvs(iter,:) = rvr;
        AY = AY * tU; SX = SX * tU;
    end
    % display iter info
    xtrm(iter) = Lm/k;
    chgv(iter) = chg_rvr;
    if idisp; fprintf(disp_str,iter,Lm/k,chg_rvr); end
    
    %% check termiating criterion
    if chg_rvr < rtol
        % kkt = AY(:,end-r+1:end) - SX(:,end-r+1:end).*(ones(m,1)*rvr');
        kkt = AY(:,end-r+1:end) - SX(:,end-r+1:end)*diag(rvr);
        kktcheck = sqrt((kkt.^2)'*ones(m,1));
        kktcheck = max(kktcheck)/max(tol,rvr(end));
        if kktcheck < ptol;  break;  end
        kktc(iter) = kktcheck;
    else
        if iter == 1; kktc(iter) = inf; else kktc(iter) = kktc(iter-1); end
    end
    
    %% look-back optimization
    xtrm(iter) = Lm / k;
    if Lm == 0; continue; end
    Xm(:,1:k) = X;
    Ym(:,1:k) = Y;
    % projection
    Im = k+1:k+Lm;
    T = X'*Xm(:,Im);
    Px = Xm(:,Im) - X*T;
    Py = Ym(:,Im) - Y*T;
    T = Px'*Px;
    % remove small vectors
    if Lm > 50
        dT = diag(T);
        [sdT,idx] = sort(dT,'descend');
        %csum = cumsum(sdT);
        %portion = csum/csum(end);
        %L = find(portion > .999,1)
        L = sum(sdT > 5e-8,1);
        if L < .95*Lm; %disp([L Lm])
            Lm = L;
            Icut = idx(1:Lm);
            Py = Py(:,Icut);
            T = T(Icut,Icut);
        end
    end
    % orthonormalize Px
    [U,D] = eig(T);
    ev = diag(D);
    [~,idx] = sort(ev,'ascend');
    e_tol = min(sqrt(eps),tol);
    cut = find(ev(idx) > e_tol,1);
    if isempty(cut); Lm = 0; continue; end
    Icut = idx(cut:end);
    L = Lm - cut + 1;
    dv = 1./sqrt(ev(idx(Icut)));
    T = U(:,Icut)*sparse(1:L,1:L,dv);
    % subspace optimization
    Yo = [Y, Py*T];
    T = Yo'*Yo;
    if issparse(T); T = full(T); end
    [U,D] = eig(T);
    [rv_sort,idx] = sort(diag(D),'ascend');
    Y = Yo*U(:,idx(end-k+1:end));
    Lm = max(0,round(L/k))*k;
    if iter < memo; Lm = Lm + k; end
    if Lm > 0
        Xm(:,(1:Lm)+k) = Xm(:,1:Lm);
        Ym(:,(1:Lm)+k) = Ym(:,1:Lm);
    end
    rvr0 = rvr; rvr = rv_sort(end-r+1:end);
    chg_rvr = norm(rvr-rvr0)/norm(rvr);
    hrvs(iter,:) = rvr;
end %iter
Out.X = X;
Out.Y = Y;
Out.memo = memo;
Out.iter = iter;
Out.chgv = chgv(1:iter);
Out.kktc = kktc(1:iter);
Out.xtrm = xtrm(1:iter);
Out.hrvs = hrvs(1:iter,:);

    %% Nested function
    %%    -- setting tolerance for terminating criterion
    function set_tolerance
        tmp = min(mn/40/k,1);
        qtol = eps^tmp;
        rtol = 5*max(sqrt(tol*qtol),5*eps);
        ptol = 5*max(tol,sqrt(eps));
    end

end %solver
