function [W,H,iter,elapse,HIS]=NeNMF(V,r,varargin)

% Non-negative Matrix Factorization via Nesterov's Optimal Gradient Method.
% NeNMF: Matlab Code for Efficient NMF Solver

% Reference
%  N. Guan, D. Tao, Z. Luo, and B. Yuan, "NeNMF: An Optimal Gradient Method
%  for Non-negative Matrix Factorization", IEEE Transactions on Signal
%  Processing, Vol. 60, No. 6, PP. 2882-2898, Jun. 2012. (DOI:
%  10.1109/TSP.2012.2190406)

% The model is V \approx WH, where V, W, and H are defined as follows:
% V (m x n-dimension): data matrix including n samples in m-dimension
% space;
% W (m x r-dimension): basis matrix including r bases in m-dimension space;
% H (r x n-dimension): coefficients matrix includeing n encodings in
% r-dimension space.

% Written by Naiyang Guan (ny.guan@gmail.com)
% Copyright 2010-2012 by Naiyang Guan and Dacheng Tao
% Modified at Sept. 25 2011
% Modified at Nov. 4 2012

% <Inputs>
%        V : Input data matrix (m x n)
%        r : Target low-rank
%
%        (Below are optional arguments: can be set by providing name-value pairs)
%        MAX_ITER : Maximum number of iterations. Default is 1,000.
%        MIN_ITER : Minimum number of iterations. Default is 10.
%        MAX_TIME : Maximum amount of time in seconds. Default is 100,000.
%        W_INIT : (m x r) initial value for W.
%        H_INIT : (r x n) initial value for H.
%        TYPE : Algorithm type (Default is 'PLAIN'),
%               'PLAIN' - NeNMF (min{.5*||V-W*H||_F^2,s.t.,W >= 0 and H >= 0}.),
%               'L1R' - L1-norm regularized NeNMF (min{.5*||V-W*H||_F^2+beta*||H||_1,s.t.,W >= 0 and H >= 0}.),
%               'L2R' - L2-norm regularized NeNMF (min{.5*||V-W*H||_F^2+.5*beta*||H||_F^2,s.t.,W >= 0 and H >= 0}.),
%               'MR' - manifold regularized NeNMF (min{.5*||V-W*H||_F^2+.5*beta*TR(H*Lp*H^T),s.t.,W >= 0 and H >= 0}.).
%        BETA : Tradeoff parameter over regularization term. Default is 1e-3.
%        S_MTX : Similarity matrix constructed by 'constructW'.
%        TOL : Stopping tolerance. Default is 1e-5. If you want to obtain a more accurate solution, decrease TOL and increase MAX_ITER at the same time.
%        VERBOSE : 0 (default) - No debugging information is collected.
%                  1 (debugging purpose) - History of computation is returned by 'HIS' variable.
%                  2 (debugging purpose) - History of computation is additionally printed on screen.
% <Outputs>
%        W : Obtained basis matrix (m x r).
%        H : Obtained coefficients matrix (r x n).
%        iter : Number of iterations.
%        elapse : CPU time in seconds.
%        HIS : (debugging purpose) History of computation,
%               niter - total iteration number spent for Nesterov's optimal
%               gradient method,
%               cpus - CPU seconds at iteration rounds,
%               objf - objective function values at iteration rounds,
%               prjg - projected gradient norm at iteration rounds.
%
% <Usage Examples>
%        >>V=rand(1000,500);
%        >>NeNMF(V,10);
%        >>NeNMF(V,20,'verbose',1);
%        >>NeNMF(V,30,'verbose',2,'w_init',rand(m,r));
%        >>NeNMF(V,5,'verbose',2,'tol',1e-5);

% Note: another file 'GetStopCriterion.m' should be included under the same
% directory as this code.

if ~exist('V','var'),    error('please input the sample matrix.\n');    end
if ~exist('r','var'),    error('please input the low rank.\n'); end

[m,n]=size(V);

% Default setting
MaxIter=1000;
MinIter=10;
MaxTime=100000;
W0=rand(m,r);
H0=rand(r,n);
type='PLAIN';
beta=1e-3;
tol=1e-5;
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
            case 'TYPE',        type=upper(varargin{i+1});
            case 'BETA',        beta=varargin{i+1};
            case 'S_MTX',       S=varargin{i+1};
            case 'TOL',         tol=varargin{i+1};
            case 'VERBOSE',     verbose=varargin{i+1};
            otherwise
                error(['Unrecognized option: ',varargin{i}]);
        end
    end
end

ITER_MAX=1000;      % maximum inner iteration number (Default)
ITER_MIN=10;        % minimum inner iteration number (Default)
global STOP_RULE;
STOP_RULE = 1;      % '1' for Projected gradient norm (Default)
                    % '2' for Normalized projected gradient norm
                    % '3' for Normalized KKT residual

% Initialization
W=W0; H=H0;
HVt=H*V'; HHt=H*H';
WtV=W'*V; WtW=W'*W;
GradW=W*HHt-HVt';
switch (type),
    case 'PLAIN',
        GradH=WtW*H-WtV;
    case 'L1R',
        GradH=WtW*H-WtV+beta;
    case 'L2R',
        GradH=WtW*H-WtV+beta*H;
    case 'MR',
        if ~exist('S','var'),
            error('Similarity matrix must be collected for NeNMF-MR.\n');
        end
        D=diag(sum(S)); Lp=D-S;
        if ~issparse(Lp),
            LpC=norm(Lp);   % Lipschitz constant of MR term
        else
            LpC=norm(full(Lp));
        end
        GradH=WtW*H-WtV+beta*H*Lp;
    otherwise
        error('No such algorithm (%s).\n',type);
end
init_delta=GetStopCriterion(STOP_RULE,[W',H],[GradW',GradH]);
tolH=max(tol,1e-3)*init_delta;
tolW=tolH;               % Stopping tolerance
constV=sum(sum(V.^2));
% Historical information
HIS.niter=0;
HIS.cpus=0;
switch (type),
    case 'PLAIN',
        HIS.objf=sum(sum(WtW.*HHt))-2*sum(sum(WtV.*H));
    case 'L1R',
        HIS.objf=sum(sum(WtW.*HHt))-2*sum(sum(WtV.*H))+2*beta*sum(sum(H));
    case 'L2R',
        HIS.objf=sum(sum(WtW.*HHt))-2*sum(sum(WtV.*H))+beta*sum(sum(H.^2));
    case 'MR',
        HIS.objf=sum(sum(WtW.*HHt))-2*sum(sum(WtV.*H))+beta*sum(sum((H'*H).*Lp));
    otherwise,
        error('No such algorithm (%s).\n',type);
end
HIS.prjg=init_delta;

% Iterative updating
elapse=cputime;
W=W';
for iter=1:MaxIter,
    % Optimize H with W fixed
    switch (type),
        case 'PLAIN',
            [H,iterH]=NNLS(H,WtW,WtV,ITER_MIN,ITER_MAX,tolH);
        case 'L1R',
            [H,iterH]=NNLS_L1(H,WtW,WtV,beta,ITER_MIN,ITER_MAX,tolH);
        case 'L2R',
            [H,iterH]=NNLS_L2(H,WtW,WtV,beta,ITER_MIN,ITER_MAX,tolH);
        case 'MR',
            [H,iterH]=NNLS_MR(H,WtW,WtV,beta,Lp,LpC,ITER_MIN,ITER_MAX,tolH);
        otherwise,
            error('No such algorithm (%s).\n',type);
    end
    if iterH<=ITER_MIN,
        tolH=tolH/10;
    end
    HHt=H*H';   HVt=H*V';
    
    % Optimize W with H fixed
    [W,iterW,GradW]=NNLS(W,HHt,HVt,ITER_MIN,ITER_MAX,tolW);
    if iterW<=ITER_MIN,
        tolW=tolW/10;
    end
    WtW=W*W'; WtV=W*V;
    switch (type),
        case 'PLAIN',
            GradH=WtW*H-WtV;
        case 'L1R',
            GradH=WtW*H-WtV+beta;
        case 'L2R',
            GradH=WtW*H-WtV+beta*H;
        case 'MR',
            GradH=WtW*H-WtV+beta*H*Lp;
        otherwise
            error('No such algorithm (%s).\n',type);
    end
    HIS.niter=HIS.niter+iterH+iterW;
    delta=GetStopCriterion(STOP_RULE,[W,H],[GradW,GradH]);
    
    % Output running detials
    if verbose,
        switch (type),
            case 'PLAIN',
                objf=sum(sum(WtW.*HHt))-2*sum(sum(WtV.*H));
            case 'L1R',
                objf=sum(sum(WtW.*HHt))-2*sum(sum(WtV.*H))+2*beta*sum(sum(H));
            case 'L2R',
                objf=sum(sum(WtW.*HHt))-2*sum(sum(WtV.*H))+beta*sum(sum(H.^2));
            case 'MR',
                objf=sum(sum(WtW.*HHt))-2*sum(sum(WtV.*H))+beta*sum(sum((H'*H).*Lp));
            otherwise,
                error('No such algorithm (%s).\n',type);
        end
        HIS.objf=[HIS.objf,objf];
        HIS.cpus=[HIS.cpus,cputime-elapse];
        HIS.prjg=[HIS.prjg,delta];
        if (verbose==2) && (rem(iter,10)==0),
            fprintf('%d:\tstopping criteria = %e,\tobjective value = %f.\n', iter,delta/init_delta,HIS.objf(end)+constV);
        end
    end
    
    % Stopping condition
    if (delta<=tol*init_delta && iter>=MinIter) || HIS.cpus(end)>=MaxTime,
        break;
    end
end
W=W';
elapse=cputime-elapse;

if verbose,
    HIS.objf=0.5*(HIS.objf+constV);
    if verbose==2,
        fprintf('\nFinal Iter = %d,\tFinal Elapse = %f.\n', iter,elapse);
    end
end

return;

% Non-negative Least Squares with Nesterov's Optimal Gradient Method
function [H,iter,Grad]=NNLS(Z,WtW,WtV,iterMin,iterMax,tol)

global STOP_RULE;

if ~issparse(WtW),
    L=norm(WtW);	% Lipschitz constant
else
    L=norm(full(WtW));
end
H=Z;    % Initialization
Grad=WtW*Z-WtV;     % Gradient
alpha1=1;

for iter=1:iterMax,
    H0=H;
    H=max(Z-Grad/L,0);    % Calculate sequence 'Y'
    alpha2=0.5*(1+sqrt(1+4*alpha1^2));
    Z=H+((alpha1-1)/alpha2)*(H-H0);
    alpha1=alpha2;
    Grad=WtW*Z-WtV;
    
    % Stopping criteria
    if iter>=iterMin,
        % Lin's stopping condition
        pgn=GetStopCriterion(STOP_RULE,Z,Grad);
        if pgn<=tol,
            break;
        end
    end
end

Grad=WtW*H-WtV;

return;

% L1-norm Regularized Non-negative Least Squares with Nesterov's Optimal Gradient Method
function [H,iter,Grad]=NNLS_L1(Z,WtW,WtV,beta,iterMin,iterMax,tol)

global STOP_RULE;

if ~issparse(WtW),
    L=norm(WtW);	% Lipschitz constant
else
    L=norm(full(WtW));
end
H=Z;    % Initialization
Grad=WtW*Z-WtV+beta;     % Gradient
alpha1=1;

for iter=1:iterMax,
    H0=H;
    H=max(Z-Grad/L,0);    % Calculate sequence 'Y'
    alpha2=0.5*(1+sqrt(1+4*alpha1^2));
    Z=H+((alpha1-1)/alpha2)*(H-H0);
    alpha1=alpha2;
    Grad=WtW*Z-WtV+beta;
    
    % Stopping criteria
    if iter>=iterMin,
        pgn=GetStopCriterion(STOP_RULE,Z,Grad);
        if pgn<=tol,
            break;
        end
    end
end

Grad=WtW*H-WtV+beta;

return;

% L2-norm Regularized Non-negative Least Squares with Nesterov's Optimal Gradient Method
function [H,iter,Grad]=NNLS_L2(Z,WtW,WtV,beta,iterMin,iterMax,tol)

global STOP_RULE;

if ~issparse(WtW),
    L=norm(WtW)+beta;	% Lipschitz constant
else
    L=norm(full(WtW))+beta;
end
H=Z;    % Initialization
Grad=WtW*Z-WtV+beta*Z;     % Gradient
alpha1=1;

for iter=1:iterMax,
    H0=H;
    H=max(Z-Grad/L,0);    % Calculate sequence 'Y'
    alpha2=0.5*(1+sqrt(1+4*alpha1^2));
    Z=H+((alpha1-1)/alpha2)*(H-H0);
    alpha1=alpha2;
    Grad=WtW*Z-WtV+beta*Z;
    
    % Stopping criteria
    if iter>=iterMin,
        pgn=GetStopCriterion(STOP_RULE,Z,Grad);
        if pgn<=tol,
            break;
        end
    end
end

Grad=WtW*H-WtV+beta*H;

return;

% Manifold Regularized Non-negative Least Squares with Nesterov's Optimal Gradient Method
function [H,iter,Grad]=NNLS_MR(Z,WtW,WtV,beta,Lp,LpC,iterMin,iterMax,tol)

global STOP_RULE;

if ~issparse(WtW),
    L=norm(WtW)+beta*LpC;	% Lipschitz constant
else
    L=norm(full(WtW))+beta*LpC;
end
H=Z;    % Initialization
Grad=WtW*Z-WtV+beta*Z*Lp;     % Gradient
alpha1=1;

for iter=1:iterMax,
    H0=H;
    H=max(Z-Grad/L,0);    % Calculate sequence 'Y'
    alpha2=0.5*(1+sqrt(1+4*alpha1^2));
    Z=H+((alpha1-1)/alpha2)*(H-H0);
    alpha1=alpha2;
    Grad=WtW*Z-WtV+beta*Z*Lp;
    
    % Stopping criteria
    if iter>=iterMin,
        pgn=GetStopCriterion(STOP_RULE,Z,Grad);
        if pgn<=tol,
            break;
        end
    end
end

Grad=WtW*H-WtV+beta*H*Lp;

return;