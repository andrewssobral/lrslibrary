%% Demo: synthetic data. 

clear all
close all hidden

addpath('./Riem_RPCP');
load_sig=0;                     % loading signal
input_sig=1;                    % 1: synthetic data
                                % 2: surveillance video

n1=5e2; n2=5e2;                 % matrix size
data.r=round(.1*min(n1,n2));   % presumed matrix rank
data.s=round(.1*n1*n2);        % presumed matrix sparsity

noise=0e-4;                     % noise level
data_formation

if ~load_sig, A_ref=[]; B_ref=[]; end


%% model setting
model_formation


%% algorithm parameters
slv.kmax=100;           % AMS max iteration
slv.ktol=1e-9;          % tolerance on residual norm

ls.imax=20;             % path search max iteration
% ls.c=.001;

cg.imax=30;             % CG max iteration
cg.tol=1e-2;            % CG residual tolerance

trim.sig=2;             % trimming signal (=2 recommended)
trim.disp=0;
trim.tol=2e-1;
trim.mod=1;             % frequency of trimming
trim.thr1=log10(5);     % threshold for singular values of A
trim.thr2=A_max*.2;     % threshold for absolute values of B

sg.sig=1;               % safeguard signal:
                        % 0: no safteguard
                        % 1: safeguard for B-subprob only (A-subprob via projected dogleg)
                        % 2: safeguard for both subproblems
sg.thr=1e-1;            % threshold for safeguard

riop.sig=2;             % 1: Riemannian gradient
                        % 2: projected dogleg
riop.eps_c=0e-6;


%% initialization
n3=round(data.r*1.5);   % initial rank of A
n4=round(data.s*1.5);   % initial cardinality of B
tic;
initialization
time_main=toc;

% return


%% alternating minimization scheme
slv_lrs_ams


%% print results
plot_stat
fprintf('final relative error=%1.5e\n',norm(A(:)-A_0(:))/norm(A_0(:)))

% save([datestr(now,'yymmddHHMMSSFFF'),'.mat'],'Z','A','B','U','S','V','stat',...
%     'slv','trim','cg','trim','sg','A_max','A_0','B_0','A_ref','B_ref')
save matlab A_0 B_0 Z A_max A_ref B_ref

