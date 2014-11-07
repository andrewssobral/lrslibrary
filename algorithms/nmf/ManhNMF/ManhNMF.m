function [W,H,iter,elapse,HIS]=ManhNMF(X,r,varargin)

% Manhattan Non-negative Matrix Factorization.
% ManhNMF: Matlab Code for Efficient Robust Manhattan NMF Solver

% Reference
%  [1] N. Guan, D. Tao, Z. Luo, and J. Shawe-taylor, "MahNMF: Manhattan
%  Non-negative Matrix Factorization," arXiv:1207.3438v1, 2012.
%  [2] N. Guan, D. Tao, Z. Luo, and J. Shawe-taylor, "MahNMF: Manhattan
%  Non-negative Matrix Factorization," Submitted to Journal of Machine Learning Research, 2013.

% The model is X \approx W^TH, where X, W, and H are defined as follows:
% X (m x n): data matrix including n samples in m-dimensional space;
% W (r x m): basis matrix including r bases in m-dimensional space;
% H (r x n): coefficients matrix includeing n encodings in r-dimensional space.

% Written by Naiyang Guan (ny.guan@gmail.com)
% Copyright 2012-2014 by Naiyang Guan and Dacheng Tao
% Modified at Jan. 28 2013

% <Inputs>
%        X : Input data matrix (m x n)
%        r : Target low-rank
%
%        (Below are optional arguments: can be set by providing name-value pairs)
%        MAX_ITER : Maximum number of iterations. Default is 1,000.
%        MIN_ITER : Minimum number of iterations. Default is 10.
%        MAX_TIME : Maximum amount of time in seconds. Default is 100,000.
%        W_INIT : (m x r) initial value for W.
%        H_INIT : (r x n) initial value for H.
%        LAM_INIT : initial value of smoothness parameter. Default is 1.
%        MDL_TYPE : Model type (Default is 'PLAIN'),
%               'PLAIN' - MahNMF (min{||X-W^T*H||_1,s.t.,W >= 0 and H >= 0}.),
%               'BXC' - Box Constrained MahNMF (min{||X-W^T*H||_1,s.t.,1 >= W >= 0 and 1 >= H >= 0}.),
%               'MNR' - Manifold Regularized MahNMF
%               (min{||X-W^T*H||_1+.5*beta*TR(H*Lp*H^T),s.t.,W >= 0 and H >= 0}.),
%               'GSP' - Group Sparse MahNMF
%               (min{||X-W^T*H||_1+.5*beta*\sum_{g\in G}||W^[g]||_{1,p},s.t.,W >= 0 and H >= 0}.),
%               'SYM' - Symmetric MahNMF (min{||X-H*H^T||_1,s.t., H >= 0}.).
%        ALG_TYPE : Algorithm type (Default is 'AGD'),
%               'AGD' - Accelerated Gradient Descent,
%               'RRI' - Rank-one Residue Iteration.
%        BETA : Tradeoff parameter over regularization term. Default is 1e-3.
%        SIM_MTX : Similarity matrix constructed by 'constructW'.
%        GPP_MTX : Group pattern for boundary of all groups.
%        TOL_INNR : Stopping tolerance of inner iterations. Default is 1e-2.
%        TOL_OUTR : Stopping tolerance of outer iterations. Default is 1e-3.
%               If you want to obtain a more accurate solution, decrease TOL_INNR or TOL_OUTR and increase MAX_ITER at the same time.
%        VB_OUTR : 0 (default) - No debugging information is collected.
%                  1 (debugging purpose) - History of computation is returned by 'HIS' variable.
%                  2 (debugging purpose) - History of computation is additionally printed on screen.
%        VB_INNR : 0 (default) - No debugging information is collected.
%                  1 (debugging purpose) - History of computation is returned by 'HIS' variable.
%                  2 (debugging purpose) - History of computation is additionally printed on screen.
% <Outputs>
%        W : Obtained basis matrix (r x m).
%        H : Obtained coefficients matrix (r x n).
%        iter : Number of iterations.
%        elapse : CPU time in seconds.
%        HIS : (debugging purpose) History of computation,
%               niter - total iteration number spent for Nesterov's optimal
%               gradient method,
%               cpus - CPU seconds at iteration rounds,
%               objf - objective function values at iteration rounds,
%               dlta - stopping criteria of block coordinate descent.
%
% <Usage Examples>
%        >>X=rand(1000,500);
%        >>ManhNMF(X,10);
%        >>ManhNMF(X,20,'verbose',1);
%        >>ManhNMF(X,30,'verbose',2,'w_init',rand(r,m));
%        >>ManhNMF(X,5,'verbose',2,'tol_outr',1e-5);

% Note: other files 'GetStopCriterion.m', 'ApproxFunC.m', and 'wmedianf.mexw32' should be included under the same
% directory as this code.

if ~exist('X','var'),    error('please input the sample matrix.\n');    end
if ~exist('r','var'),    error('please input the low rank.\n'); end

[m,n]=size(X);

% Default setting
MaxIter=1000;
MinIter=10;
MaxTime=10000;
W0=rand(r,m);
H0=rand(r,n);
lambda0 = 1;
mdl_type='PLAIN';
alg_type='AGD';
beta=1e-3;
tol_outr=1e-3;
tol_innr=1e-2;
vb_outr=0;
vb_innr=0;

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
            case 'LAM_INIT',    lambda0=varargin{i+1};
            case 'MDL_TYPE',    mdl_type=upper(varargin{i+1});
            case 'ALG_TYPE',    alg_type=upper(varargin{i+1});
            case 'BETA',        beta=varargin{i+1};
            case 'SIM_MTX',     S=varargin{i+1};
            case 'TOL_OUTR',    tol_outr=varargin{i+1};
            case 'TOL_INNR',    tol_innr=varargin{i+1};
            case 'VB_OUTR',     vb_outr=varargin{i+1};
            case 'VB_INNR',     vb_innr=varargin{i+1};
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
switch (alg_type),
    case 'AGD',
        W=W0; H=H0;
    case 'RRI',
        W=W0; H=H0;
        fprintf('NRRI: please make sure that the initial points contain small entries.\n');
    otherwise,
        error('No such algorithm (%s).\n',alg_type);
end
switch (mdl_type),
    case 'PLAIN',
    case 'BXC',
    case 'MNR',
        if ~exist('S','var'),
            error('Similarity matrix must be collected for NeNMF-MR.\n');
        end
        D=diag(sum(S)); Lp=D-S;
        if ~issparse(Lp),
            LpC=norm(Lp);   % Lipschitz constant of MR term
        else
            LpC=norm(full(Lp));
        end
    case 'GSP',
    case 'SYM',
    otherwise
        error('No such model (%s).\n',mdl_type);
end
% Historical information
HIS.objf=0;         % Historical objective values
HIS.niter=0;        % Historical iteration numbers
HIS.cpus=0;         % Historical CPU seconds
switch (mdl_type),
    case 'PLAIN',
        E = X-W'*H;          % Residue matrix
        HIS.objf(1)=sum(abs(E(:)));
    case 'BXC',
        E = X-W'*H;
        HIS.objf(1)=sum(abs(E(:)));
    case 'MNR',
        E = X-W'*H;
        HIS.objf(1)=sum(abs(E(:)))+.5*beta*sum(sum((H'*H).*Lp));
    case 'GSP',
        E = X-W'*H;
        HIS.objf(1)=sum(abs(E(:)));
    case 'SYM',
        H = H';             % W=H' in this case
        E = X-H*H';
        HIS.objf(1)=sum(abs(E(:)));
    otherwise,
        error('No such model (%s).\n',mdl_type);
end
if vb_outr==2,
    fprintf('\nManhNMF (%s) by %s starts,\tinitial objective value = %f.\n', mdl_type,alg_type,HIS.objf(1));
end

% Iterative updating
tolH=tol_innr;
tolW=tolH;
elapse=cputime;
for iter=1:MaxIter,
    % Optimize H with W fixed
    if vb_outr==2,
        fprintf('\tOptimize H ...\n');
    end
    switch (mdl_type),
        case 'PLAIN',
            switch (alg_type),
                case 'AGD',
                    [H,E,iterH]=NLAD_AGD(X,W',H,lambda0/sqrt(iter),ITER_MIN,ITER_MAX,tolH,vb_innr);
                case 'RRI',
                    [H,E,iterH]=NRRI_OPS(W',H,E,ITER_MIN,ITER_MAX,tolH,vb_innr);
                otherwise,
                    error('No such algorithm (%s).\n',alg_type);
            end
        case 'BXC',
            fprintf('The box-constrained ManhNMF code is coming.\n');
        case 'MNR',
            fprintf('The manifold regularized ManhNMF code is coming.\n');
        case 'GSP',
            fprintf('The group sparse ManhNMF code is coming.\n');
        case 'SYM',
            fprintf('The symmetric ManhNMF code is coming.\n');
        otherwise,
            error('No such model (%s).\n',mdl_type);
    end
    if iterH<=ITER_MIN,
        tolH=tolH/10;
    end
    
    % Optimize W with H fixed
    if vb_outr==2,
        fprintf('\tOptimize W ...\n');
    end
    switch (mdl_type),
        case 'PLAIN',
            switch (alg_type),
                case 'AGD',
                    [W,E,iterW]=NLAD_AGD(X',H',W,lambda0/sqrt(iter),ITER_MIN,ITER_MAX,tolW,vb_innr);   E=E';
                case 'RRI',
                    [W,E,iterW]=NRRI_OPS(H',W,E',ITER_MIN,ITER_MAX,tolW,vb_innr); E=E';
                otherwise,
                    error('No such algorithm (%s).\n',alg_type);
            end
        case 'BXC',
            fprintf('The box-constrained ManhNMF code is coming.\n');
        case 'MNR',
            fprintf('The manifold regularized ManhNMF code is coming.\n');
        case 'GSP',
            fprintf('The group sparse ManhNMF code is coming.\n');
        case 'SYM',
            fprintf('The symmetri ManhNMF code is coming.\n');
        otherwise
            error('No such model (%s).\n',mdl_type);
    end
    if iterW<=ITER_MIN,
        tolW=tolW/10;
    end
    HIS.niter=HIS.niter+iterH+iterW;
    
    % Output running detials
    switch (mdl_type),
        case 'PLAIN',
            objf=sum(abs(E(:)));
        case 'BXC',
            objf=sum(abs(E(:)));
        case 'MNR',
            objf=sum(abs(E(:)))+.5*beta*sum(sum((H'*H).*Lp));
        case 'GSP',
            objf=sum(abs(E(:)));
        case 'SYM',
            H=H';             % W=H' in this case
            E=X-H*H';
            objf=sum(abs(E(:)));
        otherwise,
            error('No such model (%s).\n',mdl_type);
    end
    HIS.objf=[HIS.objf,objf];
    HIS.cpus=[HIS.cpus,cputime-elapse];
    delta=abs(HIS.objf(iter)-HIS.objf(iter+1))/abs(HIS.objf(1)-HIS.objf(iter+1));
    HIS.dlta=delta;
    if vb_outr,
        if (vb_outr==2) && (rem(iter,1)==0),
            fprintf('%d:\titerH = %d,\titerW = %d,\tstopping criteria = %e,\tobjective value = %f.\n', iter,iterH,iterW,delta,HIS.objf(iter+1));
        end
    end
    
    % Stopping condition
    if (delta<=tol_outr && iter>=MinIter) || (HIS.cpus(end)>=MaxTime),
        break;
    end
end
elapse=cputime-elapse;

if vb_outr==2,
    fprintf('\nFinal Iter = %d,\tFinal Elapse = %f.\n', iter,elapse);
end

return;

% Non-negative Least Absolute Deviations with Nesterov's Smoothing+OGM
% Model: min{||V-WH||_1, s.t., H >= 0} with W >= 0 fixed.
function [H,E,iter]=NLAD_AGD(V,W,H,lambda,iterMin,iterMax,tol,verbose)

global STOP_RULE;

%--------- Algorithm Setting ---------
[~,n]=size(H);
d=max(eps,sqrt(sum(W.^2,2)));
D=d*ones(1,n);
L=sum(d)/lambda;                      % Lipschitz constant
precision=min(.1,tol);

%------------ OGM Setting ------------
E=V-W*H;                              % Residue at H0
U=min(1,max(-1,-E./D/lambda));        % Dual variable at H0
Grad=W'*U;                            % Gradient at H0
init_grad=GetStopCriterion(STOP_RULE,H,Grad);
SumGrad=H;                            % Prox-func center: H0
k=0;
if verbose,
    objHis=ApproxFunC(E,D,lambda);    % Objective value
    if verbose==2,
        fprintf('\t[%d]:\t%f,\tdelta=%f.\n',1,objHis,1);
    end
end

for iter=1:iterMax,
    %--------- Compute Dual, Gradient and Descent ---------
    U=min(1,max(-1,-E./D/lambda));    % Dual variable at Z
    Grad=W'*U;                        % Gradient at current Z
    delta=GetStopCriterion(STOP_RULE,H,Grad)/init_grad;
    if (delta<=precision) && (iter>=iterMin),
        break;
    end
    
    %--------- Generate 'Y' Sequence and Checking ---------
    Y=max(0,H-Grad/L);                % Get good 'Y' sequence
    SumGrad=SumGrad-Grad*(k+1)/2/L;   % Sum up hisorical gradients
    Z=max(0,SumGrad);                 % Get good 'Z' sequence
    H0=H;
    H=Z*2/(k+3)+Y*(k+1)/(k+3);        % Update H by combination
    k=k+1;                            % Update Iteration Counter
    
    %----------- Adapative Restarting Strategy ----------
    AdaRSC=sum(sum(Grad.*(H-H0)));
    if AdaRSC>0,
        k=0;
        SumGrad=H0;                   % Set to background
        H=H0;                         % Backtracking
    end
    
    %---------------- Algorithmic Updating ---------------
    E=V-W*H;                          % Update residue    
    if verbose,
        objHis(iter+1)=ApproxFunC(E,D,lambda);
        if (verbose==2) && (rem(iter,100)==0),
            fprintf('\t[%d]:\t%f,\tdelta=%f.\n',iter+1,objHis(iter+1),delta);
        end
    end
end

if verbose==2,
    fprintf('\tNLAD+OGM (%d): (L=%.2e, lam=%.2e, delta=%.3f, tol=%.2e).\n', iter,L,lambda,delta,precision);
end

return;

% Non-negativity Constrained Rank-one Residue Iteration for Subproblem of ManhNMF with OPS Weighted
% Median Solver
% Model: min{||V-WH||_1, s.t., H>=0} with W>=0 fixed.
function [H,E,iter]=NRRI_OPS(W,H,E,iterMin,iterMax,tol,verbose)

[r,n]=size(H);
objHis=sum(abs(E(:)));
if verbose==2,
    fprintf('\t[%d]:\t%f,\tdelta=%f.\n',1,objHis,1);
end

for iter=1:iterMax,
    for l=1:r,
        E=E+W(:,l)*H(l,:);    % residual eliminating rank 'l'
        row_idx=W(:,l)>0;     % filtering zero weights
        q=sum(row_idx);
        if q > 0,
            w=W(row_idx,l);
            for j=1:n,
                x=E(row_idx,j)./W(row_idx,l);
                y=wmedianf(x,w);  % weighted median
                H(l,j)=max(0,y);
            end
        end        
        E=E-W(:,l)*H(l,:);    % update the residual
    end
    
    objHis(iter+1)=sum(abs(E(:)));
    delta=abs(objHis(iter)-objHis(iter+1))/abs(objHis(iter)-objHis(1));
    if (delta<=tol) && (iter>=iterMin),
        break;
    end
    
    if (verbose==2) && (rem(iter,100)==0),
        fprintf('\t[%d]:\t%f,\tdelta=%f.\n',iter+1,objHis(iter+1),delta);
    end
end

if verbose==2,
    fprintf('\tNRRI+OPS sucesses in %d iterations, tolerance=%f.\n', iter,tol);
end

return;