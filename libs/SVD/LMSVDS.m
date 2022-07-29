function [U,S,V,Out] = LMSVDS(A,r,opts)
%
% Limited Memory Block Krylov Subspace Optimization for
%    Computing Dominant Singular Value Decompositions
%
% Input:
%    A --- either a numeric matrix or a struct with 4 fields:
%        A.size  -- [m, n]
%        A.times -- a function handle for A*x
%        A.trans -- a function handle for A'*x
%        A.param -- parameters used by A.times and A.trans
%                   (assign [] if none)
%        r --- number of leading singular triplets
%   opts --- option structure with fields:
%        tol     -- tolerance, default is 1.e-8
%        maxit   -- maximal number of iteration, default is 300
%        sub     -- subspace type, chosen from [1,2], default is adaptively
%                   chosen from 1 and 2 by the relation between r and m,n
%        memo    -- number of block subspaces, default is adaptively chosen
%                   from 5 to 3 by the relation between r and m,n
%        gvk     -- number of addional guard vectors, default is 10
%        initY   -- initial guess of a n*r matrix, default is randomly
%                   generated under Gaussian distribution
%        idisp   -- detailed information display option, default is 0
% Output:
%    U,S,V --- dominant part of the SVD of A with the r singular triplets
%    Out   --- output information
%
% Copyright 2012. Version 3.0 (v3.0a)
% Xin Liu, Yin Zhang, Zaiwen Wen. Jul. 28, 2010.
% Last Revision: March 15, 2012.
%

if nargin < 2; r = 6;  end
[m,n] = check_matrix(A);
% if r > min(m,n)/2;
%     warning(TooManySVSrequested,'r > min(m,n)/2');
% end

r = uint8(r);
% disp(['r=', num2str(r)])

% set parameters
mainargin = nargin;
tol = 1e-8;
maxit = 300;
idisp = 0;
mn = min(m,n);

if r <= mn*0.02
    memo = 5;
elseif r <= mn*0.03
    memo = 4;
else
    memo = 3;
end

if r >= mn*0.015 && r <= mn*0.035
    sub = 1;
else
    sub = 2;
end

if isfield(opts,  'gvk'); tau = opts.gvk; else tau = 10; end
% working size
k = min([2*r,r+tau,m,n]);
% initial guess
% disp(n)
% disp(['k=', num2str(k)])
% Y = randn(n,k);
Y = randn(n,round(k));

if mainargin < 3; return; end

if isfield(opts,  'tol');     tol = opts.tol;   end
if isfield(opts,'maxit');   maxit = opts.maxit; end
if isfield(opts, 'memo');    memo = opts.memo;  end
if isfield(opts,'idisp');   idisp = opts.idisp; end
if isfield(opts,'initY');       Y = opts.initY(:,1:k); end
if isfield(opts,  'sub');     sub = opts.sub;   end

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
if sub == 1
    [X,Y,Out] = lm_lbo1(A,X,Y,r,tol,maxit,memo,idisp);
else
    [X,Y,Out] = lm_lbo2(A,X,Y,r,tol,maxit,memo,idisp);
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [U,S,V] = get_svd(X,Y)
        method = 2;
        switch method
            case 1
                [V,S,W] = svd(Y,0);
                U = X*W;
            case 2
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

%%%%%%%%%%%%%%%%%%%%%% solver 1 %%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,Out] = lm_lbo1(A,X,Y,r,tol,maxit,memo,idisp)
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

m = size(X,1);
n = size(Y,1);
mn = min(m,n);
k = size(Y,2);
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

%% Stopping Criterion and Tolerance
tmp = min(mn/40/k,1);
qtol = eps^tmp;
rtol = 5*max(sqrt(tol*qtol),5*eps);
ptol = 5*max(tol,sqrt(eps));
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
    
    %% check stopping
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
        if L < .95*Lm %disp([L Lm])
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

end %solver

%%%%%%%%%%%%%%%%%%%%%% solver 2 %%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,Out] = lm_lbo2(A,X,Y,r,tol,maxit,memo,idisp)
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
%       r --- number of sv's requested

m = size(X,1);
n = size(Y,1);
mn = min(m,n);
k = size(Y,2);
% disp(['k=', num2str(k), ' : r=', num2str(r)])
if k < r; error('working size too small'); end

if memo > 0
    Xm = zeros(m,memo*k);
    Ym = zeros(n,memo*k);
    Lm = k;
    Xm(:,1:k) = X;
    Ym(:,1:k) = Y;
else
    Lm = 0;
end
% disp(r)
% rvr = zeros(r,1);
rvr = zeros(round(r),1);
chg_rvr = 1;
chgv = zeros(maxit,1);
xtrm = zeros(maxit,1);
% hrvs = zeros(maxit,r);
hrvs = zeros(maxit,round(r));
kktc = zeros(maxit,1);
disp_str = 'iter %3i: memo used %i, chg_rvr %8.4e\n';

%% Stopping Criterion and Tolerance
tmp = min(mn/40/k,1);
qtol = eps^tmp;
rtol = 5*max(sqrt(tol*qtol),5*eps);
ptol = 5*max(tol,sqrt(eps));
%% ---------------------------------
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
    
    %% check stopping
    if Lm == 0 || iter <= 3
        SYTY = SX'*AY; SYTY = 0.5*(SYTY+SYTY');
        [tU,tE] = eig(SYTY);
        rvr0 = rvr;
        rv_sort = sort(diag(tE),'ascend');
        rvr = rv_sort(end-round(r)+1:end);
        chg_rvr = norm(rvr0-rvr)/norm(rvr);
        hrvs(iter,:) = rvr;
        AY = AY * tU; SX = SX * tU;
    end
    % display iter info
    xtrm(iter) = Lm/k;
    chgv(iter) = chg_rvr;
    if idisp; fprintf(disp_str,iter,Lm/k,chg_rvr); end
    
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
    if Lm == 0; continue; end
    
    % projection
    Im = 1:Lm;
    T = X'*Xm(:,Im);
    Px = Xm(:,Im) - X*T;
    Py = Ym(:,Im) - Y*T;
    T = Px'*Px;
    if Lm == k
        Xm(:,1:k) = X;
        Ym(:,1:k) = Y;
    end
    % remove small vectors
    if Lm > 50
        dT = diag(T);
        [sdT,idx] = sort(dT,'descend');
        L = sum(sdT > 5e-10,1);
        if L < .95*Lm %disp([L Lm])
            Lm = L;
            Icut = idx(1:Lm);
            Px = Px(:,Icut);
            Py = Py(:,Icut);
            T = T(Icut,Icut);
        end
    end
    % orthonormalize Px
    if issparse(T); T = full(T); end
    [U,D] = eig(T); ev = diag(D);
    [evs,idx] = sort(ev,'ascend');
    e_tol = 10 * eps;
    cut = find(evs > e_tol,1);
    if isempty(cut); Lm = 0; continue; end
    Icut = idx(cut:end);
    L = Lm - cut + 1;
    dv = 1./sqrt(evs(Icut));
    T = U(:,Icut)*sparse(1:L,1:L,dv);
    % subspace optimization
    Yo = [Y, Py*T];
    YtY = Yo'*Yo;
    if issparse(YtY); YtY = full(YtY); end
    [U,D] = eig(YtY); rv = diag(D);
    [rv_sort,idx] = sort(rv,'ascend');
    U = U(:,idx(end-k+1:end));
    Y = Yo*U; X = [X, Px*T]*U;
    Lm = max(1,round(L/k))*k;
    if iter < memo; Lm = Lm + k; end
    if Lm > k
        Xm(:,k+1:Lm) = Xm(:,1:Lm-k);
        Ym(:,k+1:Lm) = Ym(:,1:Lm-k);
        Xm(:,1:k) = X;
        Ym(:,1:k) = Y;
    end
    % r leading ritz values of AA'
    rvr0 = rvr; rvr = rv_sort(end-r+1:end);
    chg_rvr = norm(rvr-rvr0)/norm(rvr);
    hrvs(iter,:) = rvr;
end

Out.X = X;
Out.Y = Y;
Out.memo = memo;
Out.iter = iter;
Out.chgv = chgv(1:iter);
Out.kktc = kktc(1:iter);
Out.xtrm = xtrm(1:iter);
Out.hrvs = hrvs(1:iter,:);

end %solver
