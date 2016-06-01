% R2PCP: Riemannian Robust Principal Component Pursuit (Hintermüller and Wu, 2014)
% process_video('RPCA', 'R2PCP', 'dataset/demo.avi', 'output/demo_R2PCP.avi');
Z = M;
load_sig = 0; % loading signal
input_sig = 2;
noise = 0e-4; % noise level
% data_formation
if ~load_sig, A_ref=[]; B_ref=[]; end
%[z1,z2,z3] = size(Z);
z1 = params.rows;
z2 = params.cols;
z3 = size(Z,2);
n1 = z1*z2; n2 = z3;
%%% model setting
model_formation;
%%% algorithm parameters
slv.kmax=100;           % AMS max iteration
slv.ktol=1e-4;          % tolerance on residual norm
ls.imax=20;             % line search max iteration
% ls.c=.001;
cg.imax=30;             % CG max iteration
cg.tol=1e-2;            % CG residual tolerance
trim.sig=2;             % trimming signal (=2 recommended)
trim.disp=0;
trim.tol=2e-1;
trim.mod=1;             % frequency of trimming
trim.thr1=log10(10);    % threshold for singular values of A
trim.thr2=.12;          % threshold for absolute values of B
sg.sig=1;               % safeguard signal:
% 0: no safteguard
% 1: safeguard for B-subprob only (A-subprob via projected dogleg)
% 2: safeguard for both subproblems
sg.thr=1e-1;            % threshold for safeguard
riop.sig=2;             % 1: Riemannian gradient
% 2: projected dogleg
riop.eps_c=0e-6;
%%% initialization
n3 = 5;                   % initial rank of A
n4 = round(n1*n2*.1);     % initial cardinality of B
initialization;
% alternating minimization scheme
slv_lrs_ams;
L = A;
S = B;

%figure(2), imagesc(L);
%figure(3), imagesc(S);