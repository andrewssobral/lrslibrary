function test_PROPACK(PRINT)
if nargin < 1
    PRINT = true;       % whether to display vectors
end
if ~PRINT, myDisp = @(x)x; else, myDisp = @disp; end
% --- First, test the mex functions by themselves
% warning('off','PROPACK:NotUsingMex')
randn('state',1);

% --- test reorth ---
disp(' ---- testing reorth ---- ');
N = 10; M = 5; Q = randn(N,M); Q = orth(Q);
R = randn(N,1);
alpha = 0.5; % default is 0.5
method = 0;  % default is 0 (0 for MGS, 1 for GS)
warning('off','PROPACK:NotUsingMex');
[R_test,normR_new,nre] = reorth_m(Q,R,norm(R),[],alpha,method);
% test: is this an orthonormal system?
% thresh([Q,R_new/normR_new]' * [Q,R_new/normR_new],.001)

method = 0;disp('MGS');
[R_new,normR_new,nre] = reorth(Q,R,norm(R),1:5,alpha,method);
myDisp(thresh([Q,R_new/normR_new]' * [Q,R_new/normR_new],.001));
fprintf('Error: %e\n',norm( R_new - R_test ));
method = 1;disp('CGS');
[R_new,normR_new,nre] = reorth(Q,R,norm(R),1:5,alpha,method);
myDisp( thresh([Q,R_new/normR_new]' * [Q,R_new/normR_new],.001) );
fprintf('Error: %e\n',norm( R_new - R_test ));

% --- test bdsqr ---
disp(' ---- testing bdsqr ---- ');
B = randn(10);
BB = full(B); resnrm = .1;
[S,bot] = bdsqr(diag(BB),[diag(BB,-1); resnrm]);
[S1,bot1] = bdsqr_m(diag(B),[diag(B,-1); resnrm]);
myDisp( [S,S1] );
fprintf('Error: %e\n',norm( S - S1 ));



%% Not working if A is rank deficient!
% issue is in the bidiagonalization.  bdsqr works fine
% need to make opt.delta very small (and need eta <= delta )

% --- now test all of PROPACK ---
disp(' ---- testing lansvd ---- ');
%{
    Dec '09, discovered bug in PROPACK:
    if A is rank deficient, and we ask for quite a few singular values
    (or all of them), then the 2nd or 3rd singular value is spurious,
    and the corresponding singular vectors are not o.n.

    Can fix this by forcing opt.delta to be very small
%}
M = 80; N = 90;
r = 78;  % lansvd has a spurious singular value near the beginning
A = randn(M,r) * randn(r,N);  % this doesn't work
% A = randn(M,N);  % this works
opt = []; 
% opt.cgs = 0;
opt.eta = eps;
opt.delta = eps;  % delta needs to be small!
% U(:,1:4)'*U(:,1:4)

[U,S,V] = lansvd( @(x) A*x, @(y)A'*y, M, N, 51, 'L', opt );
% [U,S,V] = lansvd( @(x) A*x, @(y)A'*y, M, N, min([M,N]), 'L', opt );
% [U,S,V] = lansvd( sparse(A),70, 'L', opt );
% [U,S,V] = lanbpro( sparse(A),80,[],opt );
% [U,S,V] = svd(full(S));
[UU,SS,VV] = svd(A);  % compare with Matlab

s = diag(S); ss = diag(SS);
% K = min( length(s), 30 );
% K = M;
K = 8; Kend = length(s); K2 = Kend-8;
myDisp('PROPACK:   MATLAB:');
myDisp([ s(1:K), ss(1:K) ]);
myDisp('  ...  ');
myDisp([ s(K2:Kend), ss(K2:Kend) ]);
fprintf('discrepancy in singular values is %e\n',norm(s(1:K)-ss(1:K)));


%% internal functions
function [sigma,bnd] = bdsqr_m(alpha,beta)
k = length(alpha);
B = spdiags([alpha(:),beta(:)],[0,-1],k+1,k);
[U,S,V] = svd(full(B),0);
sigma = diag(S);
bnd = U(end,1:k)';
