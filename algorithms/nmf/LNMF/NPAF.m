function [W,H,iter,elapse,HIS]=NPAF(V,Lpos,Lneg,r,varargin)

% Non-negative Patch Alignment Framework (NPAF) with Kullback-Leibler (KL) Divergence Metric.
% NPAF: Matlab Code for An Unified NMF Framework

% Reference:
% N. Guan, D. Tao, Z. Luo, and B. Yuan, "Non-negative Patch Alignment Framework,"
% IEEE Transactions on Neural Networks, vol. 22, no. 8, pp. 1218–1230, 2011.

% The model is V \approx WH subject to constraint 'S', where V, W, and H are defined as follows:
% V (m x n-dimension): data matrix including n samples in m-dimensional space;
% Lpos (n x n-dimension): (non-negative) positive part of alignment matrix;
% Lneg (n x n-dimension): (non-negative) negative part of alignment matrix;
% W (m x r-dimension): basis matrix including r bases in m-dimensional space;
% H (r x n-dimension): coefficients matrix including n encodings in
% r-dimensional space.

% Written by Naiyang Guan (ny.guan@gmail.com)
% Copyright 2010-2012 by Naiyang Guan and Dacheng Tao

% <Inputs>
%        V : Input data matrix (m x n)
%        r : Target low-rank or reduced dimensionality
%
%        (Below are optional arguments: can be set by providing name-value pairs)
%        MAX_ITER : Maximum number of iterations. Default is 1,000.
%        MIN_ITER : Minimum number of iterations. Default is 10.
%        MAX_TIME : Maximum amount of time in seconds. Default is 100,000.
%        W_INIT : (m x r) initial value for W.
%        H_INIT : (r x n) initial value for H.
%        ALG_TYPE : Algorithm type (Default is 'MUR'),
%               'MUR' - Multiplicative update rule,
%               'SMUR' - Squqred MUR,
%               'FGD' - Fast gradient descent,
%               'MFGD' - Multiple step-size FGD.
%        BETA : Tradeoff parameter over regularization term. Default is 1e-3.
%        TOL : Stopping tolerance. Default is 1e-5. If you want to obtain a more accurate solution, decrease TOL and increase MAX_ITER at the same time.
%        VERBOSE : 0 (default) - No debugging information is collected.
%                  1 (debugging purpose) - History of computation is returned by 'HIS' variable.
%                  2 (debugging purpose) - History of computation is additionally printed on screen.
% <Outputs>
%        W : Obtained basis matrix (m x r).
%        H : Obtained coefficients matrix (r x n).
%        iter : Number of iterations performed.
%        elapse : CPU time in seconds.
%        HIS : (debugging purpose) History of computation,
%               niter - total iteration numbers spent for FGD methods,
%               cpus - CPU seconds at iteration rounds,
%               objs - objective function values at iteration rounds,
%
% <Usage Examples>
%        >>V=rand(1000,500);
%        >>options=[];
%        >>options.Metric='Euclidean';
%        >>options.NeighborMode='KNN';
%        >>options.k=5;
%        >>options.WeightMode='HeatKernel';
%        >>options.t=1;
%        >>S=constructW(V',options);
%        >>D=diag(sum(S));
%        >>NPAF(V,D,S,10);
%        >>NPAF(V,D,S,20,'verbose',1);
%        >>NPAF(V,D,S,30,'verbose',2,'w_init',rand(1000,30));
%        >>NPAF(V,D,S,5,'verbose',2,'tol',1e-5);

if ~exist('V','var'),    error('please input data matrix.\n');    end
if ~exist('Lpos','var'),    error('please input positive part of alignment matrix.\n');    end
if ~exist('Lneg','var'),    error('please input negative part of alignment matrix.\n');    end
if ~exist('r','var'),    error('please input low rank.\n'); end

[m,n]=size(V);
V=max(V,eps);   % Handle numerical problem

% Default setting
MaxIter=1000;
MinIter=10;
MaxTime=10000;
W0=rand(m,r);
H0=rand(r,n);
AlgType='MUR';
beta=1e-3;
gamma=1e-4;
tol=1e-4;
tol_fgd=1e-4;
verbose=0;

% Read optional parameters
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'MAX_ITER',    MaxIter=varargin{i+1};
            case 'MIN_ITER',    MinIter=varargin{i+1};
            case 'MAX_TIME',    MaxTime=varargin{i+1};
            case 'W_INIT',      W0=varargin{i+1};
            case 'H_INIT',      H0=varargin{i+1};
            case 'ALG_TYPE',    AlgType=upper(varargin{i+1});
            case 'BETA',        beta=varargin{i+1};
            case 'GAMMA',       gamma=varargin{i+1};
            case 'TOL',         tol=varargin{i+1};
            case 'VERBOSE',     verbose=varargin{i+1};
            otherwise
                error(['Unrecognized option: ',varargin{i}]);
        end
    end
end

% Algorithm setting
W=max(W0,eps);
H=max(H0,eps);
Lpos=beta*Lpos+gamma*eye(n);
Lneg=beta*Lneg;
L=Lpos-Lneg;    % Laplacian matrix of the aligned graph.
Z=W*H;
Z=max(Z,eps);   % Handle numerical problem

% Historical information
HIS.niter=0;
HIS.cpus=0;
HIS.objs=KLC(V,Z)+.5*sum(sum(L.*(H'*H)));

% Iterative updating
elapse=cputime;
for iter=1:MaxIter,
    % Algorithms dependings
    switch (AlgType),
        case 'MUR',     % Normal MUR
            H=H.*(H*Lneg+W'*(V./Z))./(H*Lpos+sum(W)'*ones(1,n));
            Z=W*H;
            W=W.*((V./Z)*H')./(ones(m,1)*sum(H,2)');
            Z=W*H;
        case 'SMUR',    % Squared MUR
            H=H.*sqrt((H*Lneg+W'*(V./Z))./(H*Lpos+sum(W)'*ones(1,n)));
            Z=W*H;
            W=W.*sqrt(((V./Z)*H')./(ones(m,1)*sum(H,2)'));
            Z=W*H;
        case 'FGD',     % Normal FGD
            [H,Z]=FGD_H(V,W,H,Z,Lpos,Lneg,tol_fgd);
            [W,Z]=FGD(V',H',W',Z',tol_fgd);
            W=W';   Z=Z';
        case 'MFGD',    % Multiple step-sizes MUR
            [H,Z]=MFGD_H(V,W,H,Z,Lpos,Lneg,tol_fgd);
            [W,Z]=MFGD(V',H',W',Z',tol_fgd);
            W=W';   Z=Z';
        otherwise
            error(['No such algorithm: ',AlgType]);
    end
    
    % Handle numerical problem
    W=max(W,eps);
    H=max(H,eps);
    Z=max(Z,eps);
    
    % Stopping and checking
    HIS.objs(iter+1)=KLC(V,Z)+.5*sum(sum(L.*(H'*H)));
    delta=ObjStopCriterion(HIS.objs);
    if (delta<=tol && iter>=MinIter) || HIS.cpus(end)>=MaxTime,
        break;
    end
    
    % Output details
    if (verbose==2) && (rem(iter,10)==0),
        fprintf('%d:\tstopping criteria = %e,\tobjective value = %f.\n', iter,delta,HIS.objs(end));
    end
end
elapse=cputime-elapse;

if verbose,
    if verbose==2,
        fprintf('\nFinal Iter = %d,\tFinal Elapse = %f.\n', iter,elapse);
    end
end

return;

% Get the stopping criterion value
function delta=ObjStopCriterion(objs)

if length(objs)>1,
    delta=abs((objs(end)-objs(end-1))/(objs(end)-objs(1)));
else
    delta=1;
end

return;