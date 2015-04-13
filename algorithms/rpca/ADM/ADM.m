function out = ADM(C,tau,opts)
%
% out = ADM(C,tau,opts)
% 
% This code solves the following model
%
% min_A { tau*||A(:)||_1 + ||C - A||_* }
%
% where C is the data matrix, which will be decomposed into
% a sparse matrix A and low-rank matrix C - A.
%
% tau -- a small positive parameter

% Copyright (c) Junfeng Yang, 
% Department of Mathematics, Nanjing University, Sep. 16, 2009

%% parameter setting
beta = .25/mean(abs(C(:)));
tol = 1.e-6;
maxit = 1000;
print = 0;
if isfield(opts,'beta'); beta = opts.beta; end
if isfield(opts,'tol'); tol = opts.tol; end
if isfield(opts,'maxit'); maxit = opts.maxit; end
if isfield(opts,'print'); print = opts.print; end

%% initialization
[m,n] = size(C);
A = zeros(m,n);
B = zeros(m,n);
Lambda = zeros(m,n);
if isfield(opts,'A0');  A = opts.A0; end
if isfield(opts,'B0');  B = opts.B0; end
if isfield(opts,'Lam0'); Lambda = opts.Lam0; end

%% keep record
RECORD_ERRSP = 0; RECORD_ERRLR = 0; RECORD_OBJ = 0; RECORD_RES = 0;
if isfield(opts,'Sparse'); SP = opts.Sparse; nrmSP = norm(SP,'fro'); out.errsSP = 1; RECORD_ERRSP = 1; end
if isfield(opts,'LowRank'); LR = opts.LowRank; nrmLR = norm(LR,'fro'); out.errsLR = 1; RECORD_ERRLR = 1; end
if isfield(opts,'record_obj'); RECORD_OBJ = 1; out.obj = []; end
if isfield(opts,'record_res'); RECORD_RES = 1; out.res = []; end


% main
for iter = 1:maxit
    
    nrmAB = norm([A,B],'fro');
    
    %% A - subproblem
    X = Lambda / beta + C;
    Y = X - B;
    dA = A;
    A = sign(Y) .* max(abs(Y) - tau/beta, 0);
    dA = A - dA;
    
    %% B - subprolbme
    Y = X - A;
    dB = B;
    %[U,D,VT] = mexsvd(Y,2);
    [U,D,VT] = svd(Y);
    D = diag(D);
    ind = find(D > 1/beta);
    D = diag(D(ind) - 1/beta);
    B = U(:,ind) * D * VT(ind,:);
    dB = B - dB;
    
    %% keep record
    if RECORD_ERRSP; errSP = norm(A - SP,'fro') / (1 + nrmSP); out.errsSP = [out.errsSP; errSP]; end
    if RECORD_ERRLR; errLR = norm(B - LR,'fro') / (1 + nrmLR); out.errsLR = [out.errsLR; errLR]; end
    if RECORD_OBJ;   obj = tau*norm(A(:),1) + sum(diag(D));    out.obj = [out.obj; obj];         end
    if RECORD_RES;   res = norm(A + B - C, 'fro');             out.res = [out.res; res];         end
    
    %% stopping criterion
    RelChg = norm([dA,dB],'fro') / (1 + nrmAB);
    if print, fprintf('Iter %d, RelChg %4.2e',iter,RelChg); end
    if print && RECORD_ERRSP && RECORD_ERRLR, fprintf(', errSP %4.2e, errLR %4.2e',errSP,errLR); end
    if print, fprintf('\n'); end
    if RelChg < tol, break; end
    
    %% Update Lambda
    Lambda = Lambda - beta * (A + B - C);
end

% output
out.Sparse = A;
out.LowRank = B;
out.iter = iter;
out.exit = 'Stopped by RelChg < tol';
if iter == maxit, out.exit = 'Maximum iteration reached'; end
end
