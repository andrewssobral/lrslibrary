function [A,B,loss]=cocluster(X,k,options)

% [A,B,loss]=cocluster(X,k)
%
% solves ||X-A*B'||^2 + lambda_A_k*||A||_1 + lambda_B_k*||B||_1
%
% INPUTS:
% X: input data matrix. It uses DataSet Objects (DSO) from eigenvector to
%    color according to classes (www.eigenvector.com/software, 2011)
% k: Desired number of clusters
%
%    Optional - penalties. If k is a three-vector, then k(1) is number of
%    clusters and k(2) is the actual penalty on row clusters and k(3) is
%    the penality on column clusters. Hence, this overrides the algorithmic
%    settings.
%  
%    Otherwise, the lambdas are chosen individually for row and columns AND
%    also for each cluster.
%
% Note that the default nonnegativity constrained version seems more robust
% than without nonnegativity constraints. Avoid centering if nonnegativity
% is imposed.
%
% Scaling can be done as
%
% X.data = X.data*diag(sum(X.data.^2).^(-.5));
%
% assuming a DSO.
%
% Use plotload(A,B) to plot the clusters. Use plotload(A,B,4) to color by
% the fourth class in the DSO. If you have PLS_Toolbox, use e.g. plotgui(A)
% to plot the sample clusters.
%
% Copyright, R. Bro, ND Sidiropoulos, V. Papalexakis, E. Acar, 2011
%
%

%%%%%%%%%
% OPTIONS
%%%%%%%%%
standardoptions.repetitions = 1; % How many times to refit the model to check for local minima
standardoptions.nonnegativity = 1;
standardoptions.lambda = 'auto';
if nargin<3
    options = standardoptions;
    if nargout == 1 % Just output default options
        A = standardoptions;
        return
    end
end
if isstr(options.lambda)
    if strcmpi(options.lambda,'auto'); % Select lambda automatically
    else
        warning('options.lambda not recongized')
    end
end

if length(k)==3
    [A,B,rho,loss]=cocluster2(X,k,options);
else
    for i=1:k(1)
        kk=k;
        kk(1)=1;
        [a,b,rho,loss]=cocluster2(X,kk,options);
        if i==1
            A=a;B=b;
        else
            A=[A a];
            B=[B b];
        end
        X.data = X.data-a.data*rho*b.data'; % Seems that it doesn't matter too much how one deflates
        % aa=a.data;bb=b.data;X.data = X.data-aa*pinv(aa)*X.data*(bb*pinv(bb))';
    end
end

function [A,B,rho,loss]=cocluster2(X,k,options)

% [A,B,loss]=cocluster2(X,k)
%
% solves ||X-A*B'||^2 + lambda_A*||A||_1 + lambda_B*||B||_1
%
% INPUTS:
% X: input data matrix. It uses DataSet Objects (DSO) from eigenvector to
%    color according to classes (www.eigenvector.com/software, 2011)
% k: Desired number of clusters
%
%    Optional - penalty. If k is a two-vector, then k(1) is number of clusters and
%    k(2) is magnitude increased penalty on the clusters. Hence if k(2)
%    equals 1, no change is implied. If k(2) equals 3, the penalities are
%    three times higher than they would otherwise be.
%
%    Optionel - penalties. If k is a three-vector, then k(1) is number of
%    clusters and k(2) is the actual penalty on row clusters and k(3) is
%    the penality on column clusters. Hence, this overrides the algorithmic
%    settings.
%
%
% Note that the default nonnegativity constrained version seems more robust
% than without nonnegativity constraints. Avoid centering if nonnegativity
% is imposed.
%
% Scaling can be done as
%
% X.data = X.data*diag(sum(X.data.^2).^(-.5));
%
% assuming a DSO.
%
% Use plotload(A,B) to plot the clusters. Use plotload(A,B,4) to color by
% the fourth class in the DSO. If you have PLS_Toolbox, use e.g. plotgui(A)
% to plot the sample clusters.
%
% Copyright, R. Bro, ND Sidiropoulos, V. Papalexakis, E. Acar, 2011
%
%


%%%%%%%%%
% OPTIONS
%%%%%%%%%
standardoptions.repetitions = 1; % How many times to refit the model to check for local minima
standardoptions.nonnegativity = 1;
standardoptions.lambda = 'auto';
if nargin<3
    options = standardoptions;
    if nargout == 1 % Just output default options
        A = standardoptions;
        return
    end
end
if isstr(options.lambda)
    if strcmpi(options.lambda,'auto'); % Select lambda automatically
    else
        warning('options.lambda not recongized')
    end
end
Repeats = options.repetitions ; % Repeat the model to test for local minima;
flag = options.nonnegativity; % Set to one for

showall=0; % Set to one to show fit after each iteration
MAXNUMITER  = 1000;
SMALLNUMBER = 10^-6;
[I,J] = size(X);
isitsparse=0;


if isdataset(X)
    isdso=1;
    Xold = X;
    X = X.data;
else
    isdso=0;
end


%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE METAPARAMS
%%%%%%%%%%%%%%%%%%%%%%

if length(k)<3 % Find lambda. Start at .4 of lambda_max
    lambda_Amax=maxlam(X');
    lambda_Bmax=maxlam(X);
    
    lamA=.4*lambda_Amax;
    lamB=.4*lambda_Bmax;
    F=parafac(X,k(1),[1e-7 10 0 0 NaN 0],[flag*2 flag*2]);    % Initialize with NMF (fitted using nonneg PARAFAC)
    [A,B,rho, iterations,loss,rankdef]=SMR_norm_batch(X,k(1),lamA,lamB,eps*1000,F{1},F{2});
    
    for scale = [5 1 .5 .1 .01 .001 .0001];
    go_on=1;
    it=0;
    while go_on & it < 5% Now start adjusting. First roughly (defined by the parameter scale) 
        Aold=A;Bold=B;rho_old=rho;
        lamAold=lamA;lamBold=lamB;
        if it==0;
            disp(['Finding approximate neighborhood of lambdas. Scale of search grid is ',num2str(scale)]);
        end
        
        it = it+1;
        [lamA,lamB,go_on,direction]=adjustlambda(lamA,lamB,A,B,lambda_Amax,lambda_Bmax,scale);
        if it==1;
            d=direction; 
        end% Save the first direction of change in lambda
        F=parafac(X,k(1),[1e-7 10 0 0 NaN 0],[flag*2 flag*2]);    % Initialize with NMF (fitted using nonneg PARAFAC)
        [A,B,rho, iterations,loss,rankdef]=SMR_norm_batch(X,k(1),lamA,lamB,eps*1000,F{1},F{2});
        if direction==-d % Then we are moving in opposite direction. Hence stop
            go_on=0;
            disp('OK .. stopping now. We have to refine')
        end
        if any(sum(A==0)==size(A,1))|any(sum(B==0)==size(B,1)) % step back and take the prior solution
            go_on=0;
            A=Aold;B=Bold;rho=rho_old;
            lamA=lamAold;lamB=lamBold;
        end
    end
    disp(['Penalty on clusters. Row: ',num2str(lamA),'. Column: ',num2str(lamB)])
    end
    
else
    lamA = k(2);lamB=k(3);
    F=parafac(X,k(1),[1e-7 10 0 0 NaN 0],[flag*2 flag*2]);    % Initialize with NMF (fitted using nonneg PARAFAC)
    [A,B,rho, iterations,loss,rankdef]=SMR_norm_batch(X,k(1),lamA,lamB,eps*1000,F{1},F{2});
    disp(['Penalty on clusters. Row: ',num2str(lamA),'. Column: ',num2str(lamB)])
end





if isdso % Add all the additional info from the original data for easier plotting with plotgui afterwards
    A = copydsfields(Xold,dataset(A),{1 1});
    B = copydsfields(Xold,dataset(B),{2 1});
end




function [A,B,rho, iterations,cost,rankdef]=SMR_norm_batch(X,F,lambA,lambB,SMALLNUMBER,Aold,Bold)
% Vagelis Papalexakis - 2011
if nargin==4
    SMALLNUMBER=10^-8;
end
MAXNUMITER = 10000;
quiet=1;
% [I J K]=size(X);
% A=zeros(I,F);B=zeros(J,F);
iterations = zeros(1,F);

if nargin>5
    [A B rho it cost] = SMR_norm_core(X,F,lambA,lambB,SMALLNUMBER,Aold,Bold);
else
    [A B rho it cost] = SMR_norm_core(X,F,lambA,lambB,SMALLNUMBER);
end
iterations = it;

if any(sum(abs(A))<eps*100)|any(sum(abs(B))<eps*100)
    rankdef=1;
else
    rankdef=0;
end


function [a b rho it cost] = SMR_norm_core(X,F,lambA,lambB,SMALLNUMBER,Aold,Bold)

MAXNUMITER = 10000;
if nargin>5
    a = Aold;
    b = Bold;
else
    %     [u s v] = svds(X,F);
    %     a = abs(u*s); b = abs(v);
    [I J] = size(X);
    a = rand(I,F); b = rand(J,F);
    a_s = zeros(F,1);b_s = zeros(F,1);
end
for i = 1:F
    a_s(i)=(max(a(:,i)));a(:,i)=a(:,i)/a_s(i);
    b_s(i)=(max(b(:,i)));b(:,i)=b(:,i)/b_s(i);
end

rho = diag(a_s.*b_s);
rho_bound = max(max(X));

%now begin the ALS iteration
cost = norm(X - a*rho*b' ,'fro')^2 + lambA*sum(sum(abs(a))) + lambB*sum(sum(abs(b)));
costold = 2*cost;
it = 0;
while abs((cost-costold)/costold) > SMALLNUMBER && it < MAXNUMITER && cost > 10^5*eps
    it=it+1;
    costold=cost;
    
    %re-estimate B:
    b = SMR_01(X,a*rho,b,lambB);
    %update rho
    rho = update_rho(X,a,b,rho_bound);
    % re-estimate A:
    a = SMR_01(X',b*rho,a,lambA);
    %update rho
    rho = update_rho(X,a,b,rho_bound);
    
    
    cost = norm(X - a*rho*b' ,'fro')^2 + lambA*sum(sum(abs(a))) + lambB*sum(sum(abs(b)));
    
    %         allcosts(f,it) = cost;
    
    
    %         fprintf('iteration: %d cost: %12.10f diff: %.12f\n',it,cost,abs((cost-costold)/costold));
    if costold < cost
        disp(['*** bummer! *** ',num2str(costold-cost)])
    end
end


%vec(A*diag(rho)*B')= (B krp A) * rho
%min ||vec(X)-(B krp A) * rho||^2 ---> rho = pinv(B krp A)*vec(X)
function rho = update_rho(X,A,B,rho_bound)
rho = zeros(size(A,2),size(A,2));
x = X(:);
rho_vec = pinv(krb(B,A))*x;
%rho_vec(rho_vec<0) =0;
rho_vec = min(rho_vec,rho_bound);
rho = diag(rho_vec);


function B = SMR_01(X,A,B,lambda)

[I,J]=size(X);
[I,F]=size(A);
maxit=10000;
convcrit = 1e-9;
it=0;
Oldfit=1e100;
Diff=1e100;

while Diff>convcrit && it<maxit
    it=it+1;
    for j=1:J,
        for f=1:F,
            %             data = X(:,j);
            data = X(:,j) - A*B(j,:).' + A(:,f)*B(j,f);
            alpha = A(:,f);
            
            diff= alpha'*data - lambda/2;
            if ( diff < 0)
                B(j,f) = 0;
            elseif ( diff > (alpha.'*alpha))
                B(j,f) = 1;
            else
                B(j,f)=diff/(alpha.'*alpha);
            end
            
        end
    end
    
    fit= norm(X-A*B.','fro')^2+ lambda*sum(sum(abs(B)));
    Diff=abs(Oldfit-fit);
    %     sprintf('iter %d, %.10f',it, Diff)
    Oldfit=fit;
    
end





function [A,B,loss,rankdef]=SMR_2(X,k,lambda, flag,A,B)
%solves ||X-A*B'||^2 + lambda*||A||_1 + lambda*||B||_1
%Inputs:
%X: input data matrix
%k: desired approximation rank
%lambda

showall=0; % Set to one to show fit after each iteration
MAXNUMITER  = 1000;
SMALLNUMBER = 10^-6;
[m,n] = size(X);
meth=2;

% For line search
lineparam.Delta = 4;
acc_pow=2;  % Extrapolate to the iteration^(1/acc_pow) ahead
acc_fail=0; % Indicate how many times acceleration have failed
max_fail=4; % Increase acc_pow with one after max_fail failure
lineparam.acc_pow=acc_pow;
lineparam.acc_fail=acc_fail;
lineparam.max_fail=max_fail;


new_norm=norm(X-A*B','fro')^2 + lambda*sum(sum(abs(A))) + lambda*sum(sum(abs(B)));
old_norm=new_norm+10^10;
iter_count=1;
while(abs(new_norm-old_norm)>SMALLNUMBER && iter_count < MAXNUMITER)%until convergence
    Ao=A;Bo=B;
    iter_count=iter_count+1;
    if flag
        if meth==1
            B=NNSMREW(X,A,B,lambda);
            A=NNSMREW(X',B,A,lambda);
        else
            % Speed up of Evrim
            B=NNSMREW_v1(X,A,B,lambda);
            A=NNSMREW_v1(X',B,A,lambda);
        end
    else
        if meth==1
            B = SMREW(X,A,B,lambda);
            A = SMREW(X',B,A,lambda);
        else
            B = SMREW_v1(X,A,B,lambda);
            A = SMREW_v1(X',B,A,lambda);
        end
    end
    old_norm=new_norm;
    new_norm=norm(X-A*B','fro')^2 + lambda*sum(sum(abs(A))) +lambda*sum(sum(abs(B)));
    
    if (iter_count/2 == round(iter_count/2) && iter_count>5)   % Every second iteration do a line search
        if showall
            disp('Yihaaa .... extrapolating')
        end
        lineparam.it = iter_count;
        [A,B,lineparam] = linesrch(X,A,Ao,B,Bo,lambda,lineparam);
        [new_norm] = lossval(X,A,B,lambda);
    end
    
    
    if showall|rem(iter_count,1000)==0
        %      disp(sprintf('SMR: Iteration: %d, fit: %.10f, diff: %.10f',iter_count,new_norm,abs(new_norm-old_norm)))
        lineparam.acc_pow;
    end
    
end
%disp(sprintf('SMR: Iteration: %d, fit: %.10f, diff: %.10f',iter_count,new_norm,abs(new_norm-old_norm)))
loss = new_norm;

if any(sum(abs(A))<eps*100)|any(sum(abs(B))<eps*100)
    rankdef=1;
else
    rankdef=0;
end


function [Anew,Bnew,lineparam] = linesrch(X,A,Ao,B,Bo,lambda,lineparam);

acc_pow=lineparam.acc_pow;
acc_fail=lineparam.acc_fail;
max_fail=lineparam.max_fail;
it = lineparam.it;

Fitnow = lossval(X,A,B,lambda);
acc=0;
dL = extrapol({A,B},{Ao,Bo},max(log(it),2)^(1/acc_pow));
Fitnew = lossval(X,dL{1},dL{2},lambda);

if Fitnew>Fitnow
    acc_fail=acc_fail+1;
    dL={A,B};
    if acc_fail==max_fail,
        acc_pow=acc_pow+1+1;
        acc_fail=0;
    end
end

lineparam.acc_pow=acc_pow;
lineparam.acc_fail=acc_fail;
lineparam.max_fail=max_fail;
Anew=dL{1};
Bnew=dL{2};
%---------------------------------------------------------------------
function dA = extrapol(A,Ao,delta)

dA = A;
for i=1:length(A)
    dA{i} = Ao{i}+delta*(A{i}-Ao{i});
end

function [new_norm] = lossval(X,A,B,lambda);

new_norm=norm(X-A*B','fro')^2 + lambda*sum(sum(abs(A))) + lambda*sum(sum(abs(B)));







%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%

% PARAFAC

function [Factors,it,err,corcondia]=parafac(X,Fac,Options,const,OldLoad,FixMode,Weights);

% PARAFAC multiway parafac model
%
% See also:
% 'npls' 'tucker' 'dtld' 'gram'
%
%
%     ___________________________________________________
%
%                  THE PARAFAC MODEL
%     ___________________________________________________
%
% [Factors,it,err,corcondia] = parafac(X,Fac,Options,const,OldLoad,FixMode,Weights);
%
% or skipping optional in/outputs
%
% Factors = parafac(X,Fac);
%
% Algorithm for computing an N-way PARAFAC model. Optionally
% constraints can be put on individual modes for obtaining
% orthogonal, nonnegative, or unimodal solutions. The algorithm
% also handles missing data. For details of PARAFAC
% modeling see R. Bro, Chemom. Intell. Lab. Syst., 1997.
%
% Several possibilities exist for speeding up the algorithm.
% Compressing has been incorporated, so that large arrays can be
% compressed by using Tucker (see Bro & Andersson, Chemom.
% Intell. Lab. Syst., 1998).
% Another acceleration method incorporated here is to
% extrapolate the individual loading elements a number of
% iterations ahead after a specified number of iterations.
%
% A temporary MAT-file called TEMP.mat is saved for every
% 50 iterations. IF the computer breaks down or the model
% seems to be good enough, one can break the program and
% load the last saved estimate. The loadings in TEMP.MAT
% are given a cell array as described below and can be
% converted to A, B, C etc. by FAC2LET.M typing
% [A,B,C]=fac2let(Factors,size(X));
%
% All loading vectors except in first mode are normalized,
% so that all variance is kept in the first mode (as is
% common in two-way PCA). The components are arranged as
% in PCA. After iterating, the most important component is
% made the first component etc.
%
%
%
% ----------------------INPUT---------------------
%
% X          X is the input array, which can be from three- to N-way (also
%            twoway if the third mode is interpreted as a onedimensional
%            mode).
%
% Fac        No of factors/components sought.
%
%
% ----------------OPTIONAL INPUT---------------------
%
% Options    Optional parameters. If not given or set to zero or [],
%            defaults will be used. If you want Options(5) to be 2 and
%            not change others, simply write Options(5)=2. Even if Options
%            hasn't been defined Options will contain zeros except its
%            fifth element.
%
%            Options(1) - Convergence criterion
%            The relative change in fit for which the algorithm stops.
%            Standard is 1e-6, but difficult data might require a lower value.
%
%            Options(2) - Initialization method
%            This option is ignored if PARAFAC is started with old values.
%            If no default values are given the default Options(2) is 0.
%            The advantage of using DTLD or SVD for initialization is that
%            they often provide good starting values. However, since the
%            initial values are then fixed, repeating the fitting will give
%            the exact same solution. Therefore it is not possible to substantiate
%            if a local minimum has been reached. To avoid that use an initialization
%            based on random values (2).
%
%            0  = fit using DTLD/GRAM for initialization (default if
%                                 three-way and no missing and if sizes are
%                                 largere than number of factors at least
%                                 in two modes)
%            1  = fit using SVD vectors for initialization (default if higher than three-way or missing)
%            2  = fit using random orthogonalized values for initialization
%            10 = fit using the best-fitting models of several models
%            fitted using a few iterations
%
%            Options(3) - Plotting options
%            0 = no plots
%            1 = produces several graphical outputs
%            2 = produces several graphical outputs (loadings also shown during iterations)
%            3 = as 2 but no core consistency check (very slow for large arrays and/or many components)
%
%            Options(4) - Scaling
%            0 or 1 = default scaling (columns in mode one carry the variance)
%            2      = no scaling applied (hence fixed elements will not be modified
%
%            Options(5) - How often to show fit
%            Determines how often the deviation between the model and the data
%            is shown. This is helpful for adjusting the output to the number
%            of iterations. Default is 10. If showfit is set to NaN, almost no
%            outputs are given
%
%            Options(6) - Maximal number of iterations
%            Maximal number of iterations allowed. Default is 2500.
%
% const      A vector telling type of constraints put on the loadings of the
%            different modes. Same size as DimX but the i'th element tells
%            what constraint is on that mode.
%            0 => no constraint,
%            1 => orthogonality
%            2 => nonnegativity
%            3 => unimodality (and nonnegativitiy)
%            4 => L1 fitting (will be imposed in all modes)
%            5 => L1 fitting and nonnegativity (will be imposed in all modes)
%            If const is not defined, no constraints are used.
%            For no constraints in a threeway problem const = [0 0 0]
%
% OldLoad    If initial guess of the loadings is available. OldLoad should be
%            given a cell array where OldLoad{1}=A,OldLoad{2}=B etc.
%
% FixMode    FixMode is a binary vector of same sixe as DimX. If
%            FixMode(i) = 1 => Mode i is fixed (requires old values given)
%            FixMode(i) = 0 => Mode i is not fixed hence estimated
%            Ex.: FixMode = [0 1 1] find the scores of a data set given the loadings.
%            When some modes are fixed, the numbering of the components will
%            also be fixed. Normally components are sorted according to variance
%            as in PCA, but this will not be performed if some modes are fixed.
%
% Weights    If a matrix of the same size as X is given, weighted regression
%            is performed using the weights in the matrix Weights. Statistically
%            the weights will usually contain the inverse error standard
%            deviation of the particular element
%
% ---------------------OUTPUT---------------------
%
% Factors    PARAFAC estimate of loadings in one matrix. For a 3 component
%            solution to a 4 x 3 x 3 array the loadings A, B & C will be
%            stored in a 3 element cell vector:
%            Factors{1}=A,
%            Factors{2}=B
%            Factors{3}=C
%            etc.
%
%            Use FAC2LET.M for converting to "normal" output or simply extract the
%            components as e.g. A = Factors{1};
%
% it         Number of iterations used. Can be helpful for checking if the algorithm
%            has converged or simply hit the maximal number of iterations (default 2500).
%
% err        The fit of the model = the sum of squares of errors (not including missing
%            elements).
%
% Corcondia  Core consistency test. Should ideally be 100%. If significantly below
%            100% the model is not valid
%
%
%
% OTHER STUFF
%
%  Missing values are handled by expectation maximization only. Set all
%  missing values to NaN
%
%  COMMAND LINE (SHORT)
%
%  Factors = parafac(X,Fac);
%

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.



% $ Version 1.03 $ Date 1. October   1998 $ Not compiled $ Changed sign-convention because of problems with centered data
% $ Version 1.04 $ Date 18. February 1999 $ Not compiled $ Removed auxiliary line
% $ Version 1.06 $ Date 1. December  1999 $ Not compiled $ Fixed bug in low fit error handling
% $ Version 1.07 $ Date 17. January  2000 $ Not compiled $ Fixed bug in nnls handling so that the algorithm is not stopped until nonnegative appear
% $ Version 1.08 $ Date 21. January  2000 $ Not compiled $ Changed init DTLD so that primarily negative loadings are reflected if possible
% $ Version 1.09 $ Date 30. May 2000 $ Not compiled $ changed name noptioPF to noptiopf
% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $
% $ Version 2.001 $ June 2001 $ Fixed error in weighted regression $ RB $ Not compiled $
% $ Version 2.002 $ Jan 2002 $ Fixed scaling problem due to non-identifiability of DTLD(QZ) by scaling and normalizing after each iteration $ RB $ Not compiled $
% $ Version 2.003 $ Jan 2002 $ Fixed negative solutions when nonneg imposed $ RB $ Not compiled $
% $ Version 2.004 $ Jan 2002 $ Changed initialization when many components used $ RB $ Not compiled $
% $ Version 2.005 $ Jan 2002 $ Changed absolute fit criterion (approacing eps) into relative sse/ssx$ RB $ Not compiled $
% $ Version 2.006 $ Jan 2002 $ Fixed post-scaling when fixed loadings $ RB $ Not compiled $
% $ Version 2.01 $ Jan 2003 $ Removed corcondia for two-way data (doesn't work) and fixed a bug for data with dimension 2 $ RB $ Not compiled $
% $ Version 2.011 $ feb 2003 $ Added an option (4) for not post scaling components $ RB $ Not compiled $
% $ Version 2.10  $ jan 2004 $ Fixed a plotting error occuring when fitting model to old data $ RB $ Not compiled $
% $ Version 2.11  $ jan 2004 $ Fixed that PCA can be fitted $ RB $ Not compiled $
% $ Version 2.12  $ Jul 2004 $ Fixed initialization bug $ RB $ Not compiled $
% $ Version 2.13 $ Jan 2005 $ Modified sign conventions of scores and loads $ RB $ Not compiled $
% $ Version 2.14 $ Feb 2006 $ Fixed bug in sign-swicth when loadings are fixed $ RB $ Not compiled $
% $ Version 2.15 $ Aug 2007 $ ALS scheme substantially improved using linesearch acc. to Rajih, Comon, Harshman 2007 $ RB $ Not compiled $
% $ Version 2.16 $ Jun 2010 $ Fixed that scaling of component is not done in fixed modes $ RB $ Not compiled $

NumbIteraInitia=20;
%accel_pattern='none';
%accel_pattern='shifted'; % Use every-iteration linesearch
accel_pattern='paired'; % Use every second iteration linesearch
if nargin==0
    disp(' ')
    disp(' ')
    disp(' THE PARAFAC MODEL')
    disp(' ')
    disp(' Type <<help parafac>> for more info')
    disp('  ')
    disp(' [Factors,it,err,Corcondia] = parafac(X,Fac,Options,const,OldLoad,FixMode,Weights);')
    disp(' or short')
    disp(' Factors = parafac(X,Fac);')
    disp(' ')
    disp(' Options=[Crit Init Plot NotUsed ShowFit MaxIt]')
    disp(' ')
    disp(' ')
    disp(' EXAMPLE:')
    disp(' To fit a four-component PARAFAC model to X of size 6 x 2 x 200 x 3 type')
    disp(' Factors=parafac(X,4)')
    disp(' and to obtain the scores and loadings from the output type')
    disp(' [A,B,C,D]=fac2let(Factors);')
    return
elseif nargin<2
    error(' The inputs X, and Fac must be given')
end

DimX = size(X);
X = reshape(X,DimX(1),prod(DimX(2:end)));

nonneg_obeyed = 1; % used to check if noneg is ok

if nargin<3
    OptionsDefault=[0 0 0 1 10 2500];
else
    % Call the current Options OptionsHere and load default to use if some of the current settings should be default
    Options=Options(:);
    I=length(Options);
    if I==0
        Options=zeros(8,1);
    end
    I=length(Options);
    if I<8
        Options=[Options;zeros(8-I,1)];
    end
    OptionsHere=Options;
    Options = [1e-6 0 0 1 10 2500];
    OptionsDefault=Options;
    Options=OptionsHere;
end

if ~exist('OldLoad')==1
    OldLoad=0;
elseif length(OldLoad)==0
    OldLoad=0;
end

% Convergence criteria
if Options(1,1)==0
    Options(1,1)=OptionsDefault(1,1);
end
crit=Options(1);


% Initialization
if ~any(Options(2))
    Options(2)=OptionsDefault(2);
end
Init=Options(2);

% Interim plotting
Plt=Options(3,1);
if ~any([0 1 2 3]==Plt)
    error(' Options(3,1) - Plotting - not set correct; must be 0,1,2 or 3')
end

if Options(5,1)==0
    Options(5,1)=OptionsDefault(5,1);
end
showfit=Options(5,1);
if isnan(showfit)
    showfit=-1;
end
if showfit<-1|round(showfit)~=showfit
    error(' Options(5,1) - How often to show fit - not set correct; must be positive integer or -1')
end

if Options(6,1)==0
    Options(6,1)=OptionsDefault(6);
    maxit=Options(6,1);
elseif Options(6)>0&round(Options(6))==Options(6)
    maxit=Options(6,1);
else
    error(' Options(6,1) - Maximal number of iterations - not set correct; must be positive integer')
end

ShowPhi=0; % Counter. Tuckers congruence coef/Multiple cosine/UUC shown every ShowPhiWhen'th time the fit is shown
ShowPhiWhen=10;
MissConvCrit=1e-4; % Convergence criterion for estimates of missing values
NumberOfInc=0; % Counter for indicating the number of iterations that increased the fit. ALS algorithms ALLWAYS decrease the fit, but using outside knowledge in some sense (approximate equality or iteratively reweighting might cause the algorithm to diverge

% INITIALIZE
if showfit~=-1
    disp(' ')
    disp(' PRELIMINARY')
    disp(' ')
end
ord=length(DimX);

if showfit~=-1
    disp([' A ',num2str(Fac),'-component model will be fitted'])
end

if exist('const')~=1
    const=zeros(size(DimX));
elseif length(const)~=ord
    if length(DimX)==2 & length(const)==3
        const = const(1:2);
    elseif ~iscell(const)
        const=zeros(size(DimX));
        if showfit~=-1
            disp(' Constraints are not given properly')
        end
    end
end

if showfit~=-1
    for i=1:ord
        if ~iscell(const)
            if const(i)==0
                disp([' No constraints on mode ',num2str(i)])
            elseif const(i)==1
                disp([' Orthogonality on mode ',num2str(i)])
            elseif const(i)==2
                disp([' Nonnegativity on mode ',num2str(i)])
            elseif const(i)==3
                disp([' Unimodality on mode ',num2str(i)])
            elseif any(const==4)
                disp([' L1 fitting in mode ',num2str(i)])
            elseif any(const==5)
                disp([' L1 and nonnegativity fitting in mode ',num2str(i)])
            end
        else
            disp([' compressed nonnegativity in mode ',num2str(i)])
        end
    end
end

% Check if orthogonality required on all modes
DoingPCA= 0;
if ~iscell(const)
    if max(max(const))==1
        if min(min(const))==1,
            if length(DimX)>2
                disp(' ')
                disp(' Not possible to orthogonalize all modes in this implementation.')
                error(' Contact the authors for further information')
            else
                const = [1 0]; % It's ok for PCA but do in one mode to get LS and then orthogonalize afterwards
                DoingPCA = 1;
            end
        end
    end
end

if exist('FixMode')==1
    if length(FixMode)~=ord
        FixMode = zeros(1,ord);
    end
else
    FixMode = zeros(1,ord);
end

if showfit~=-1
    if any(FixMode)
        disp([' The loadings of mode : ',num2str(find(FixMode(:)')),' are fixed'])
    end
end
if exist('Weights')~=1
    Weights=[];
end

% Display convergence criterion
if showfit~=-1
    disp([' The convergence criterion is ',num2str(crit)])
end

% Define loading as one ((r1*r2*r3*...*r7)*Fac x 1) vector [A(:);B(:);C(:);...].
% The i'th loading goes from lidx(i,1) to lidx(i,2)
lidx=[1 DimX(1)*Fac];
for i=2:ord
    lidx=[lidx;[lidx(i-1,2)+1 sum(DimX(1:i))*Fac]];
end

% Check if weighted regression required
if size(Weights,1)==size(X,1)&prod(size(Weights))/size(X,1)==size(X,2)
    Weights = reshape(Weights,size(Weights,1),prod(size(Weights))/size(X,1));
    if showfit~=-1
        disp(' Given weights will be used for weighted regression')
    end
    DoWeight=1;
else
    if showfit~=-1
        disp(' No weights given')
    end
    DoWeight=0;
end

% Make idx matrices if missing values
if any(isnan(X(:)))
    MissMeth=1;
else
    MissMeth=0;
end
if MissMeth
    id=sparse(find(isnan(X)));
    idmiss2=sparse(find(~isnan(X)));
    if showfit~=-1
        disp([' ', num2str(100*(length(id)/prod(DimX))),'% missing values']);
        disp(' Expectation maximization will be used for handling missing values')
    end
    SSX=sum(sum(X(idmiss2).^2)); % To be used for evaluating the %var explained
    % If weighting to zero should be used
    % Replace missing with mean values or model estimates initially
    %Chk format ok.
    dimisok = 1;
    if length(OldLoad)==length(DimX)
        for i=1:length(DimX)
            if ~all(size(OldLoad{i})==[DimX(i) Fac])
                dimisok = 0;
            end
        end
    else
        dimisok = 0;
    end
    if dimisok
        model=nmodel(OldLoad);
        model = reshape(model,DimX);
        X(id)=model(id);
    else
        meanX=mean(X(find(~isnan(X))));
        meanX=mean(meanX);
        X(id)=meanX*ones(size(id));
    end
else
    if showfit~=-1
        disp(' No missing values')
    end
    if iscell(const)
        SSX=sum(sum(X.^2)); % To be used for evaluating the %var explained
    elseif ~(any(const==4)|any(const==5))
        SSX=sum(sum(X.^2)); % To be used for evaluating the %var explained
    else
        SSX=sum(abs(X(:)));
    end
end

% Check if weighting is tried used together with unimodality or orthogonality
if ~iscell(const)
    if any(const==3)|any(const==1)
        if DoWeight==1
            disp(' ')
            disp(' Weighting is not possible together with unimodality and orthogonality.')
            disp(' It can be done using majorization, but has not been implemented here')
            disp(' Please contact the authors for further information')
            error
        end
    end
end

% Acceleration
acc=-5;
do_acc=1;   % Do acceleration every do_acc'th time
acc_pow=2;  % Extrapolate to the iteration^(1/acc_pow) ahead
acc_fail=0; % Indicate how many times acceleration have failed
max_fail=4; % Increase acc_pow with one after max_fail failure
if showfit~=-1
    disp(' Line-search acceleration scheme initialized')
end

% Find initial guesses for the loadings if no initial values are given

% Use old loadings
if length(OldLoad)==ord % Use old values
    if showfit~=-1
        disp(' Using old values for initialization')
    end
    Factors=OldLoad;
    % Use DTLD
elseif Init==0
    if min(DimX)>1&ord==3&MissMeth==0
        if sum(DimX<Fac)<2
            if showfit~=-1
                disp(' Using direct trilinear decomposition for initialization')
            end
            try
                [A,B,C]=dtld(reshape(X,DimX),Fac);
            catch
                A = rand(DimX(1),Fac);B = rand(DimX(2),Fac);C = rand(DimX(3),Fac);
            end
        else
            if showfit~=-1
                disp(' Using random values for initialization')
            end
            for i=1:length(DimX)
                Factors{i}=rand(DimX(i),Fac);
            end
            A = Factors{1};B=Factors{2};C = Factors{3};
        end
        A=real(A);B=real(B);C=real(C);
        % Check for signs and reflect if appropriate
        for f=1:Fac
            if sign(sum(A(:,f)))<0
                if sign(sum(B(:,f)))<0
                    B(:,f)=-B(:,f);
                    A(:,f)=-A(:,f);
                elseif sign(sum(C(:,f)))<0
                    C(:,f)=-C(:,f);
                    A(:,f)=-A(:,f);
                end
            end
            if sign(sum(B(:,f)))<0
                if sign(sum(C(:,f)))<0
                    C(:,f)=-C(:,f);
                    B(:,f)=-B(:,f);
                end
            end
        end
        Factors{1}=A;Factors{2}=B;Factors{3}=C;
        
    else
        if showfit~=-1
            disp(' Using singular values for initialization')
        end
        try
            Factors=ini(reshape(X,DimX),Fac,2);
        catch
            Factors=[];
            for i=1:length(DimX);
                l = rand(DimX(i),Fac);
                Factors{i} =l;
            end
            if showfit~=-1
                disp(' Oops sorry - ended up with random instead')
            end
            
        end
    end
    
    % Use SVD
elseif Init==1
    if all(DimX>=Fac)
        if showfit~=-1
            disp(' Using singular values for initialization')
        end
        try
            Factors=ini(reshape(X,DimX),Fac,2);
        catch
            Factors=[];
            for i=1:length(DimX);
                l = rand(DimX(i),Fac);
                Factors = [Factors;l(:)];
            end
        end
    else
        if showfit~=-1
            disp(' Using random values for initialization')
        end
        for i=1:length(DimX)
            Factors{i}=rand(DimX(i),Fac);
        end
    end
    
    % Use random (orthogonal)
elseif Init==2
    if showfit~=-1
        disp(' Using orthogonal random for initialization')
    end
    Factors=ini(reshape(X,DimX),Fac,1);
    
elseif Init==3
    error(' Initialization option set to three has been changed to 10')
    
    % Use several small ones of the above
elseif Init==10
    if showfit~=-1
        disp(' Using several small runs for initialization')
    end
    Opt=Options;
    Opt(5) = NaN;
    Opt(6) = NumbIteraInitia;
    Opt(2) = 0;
    ERR=[];
    [Factors,it,err] = parafac(reshape(X,DimX),Fac,Opt,const,[],[],Weights);
    ERR = [ERR;err];
    Opt(2) = 1;
    [F,it,Err] = parafac(reshape(X,DimX),Fac,Opt,const,[],[],Weights);
    ERR=[ERR;Err];
    if Err<err
        Factors=F;
        err=Err;
    end
    Opt(2)=2;
    for rep=1:3
        [F,it,Err]=parafac(reshape(X,DimX),Fac,Opt,const,[],[],Weights);
        ERR=[ERR;Err];
        if Err<err
            Factors=F;
            err=Err;
        end
    end
    if showfit~=-1
        disp(' ')
        disp(' Obtained fit-values')
        disp([' Method   Fit'])
        disp([' DTLD     ',num2str(ERR(1))])
        disp([' SVD      ',num2str(ERR(2))])
        disp([' RandOrth ',num2str(ERR(3))])
        disp([' RandOrth ',num2str(ERR(4))])
        disp([' RandOrth ',num2str(ERR(5))])
    end
else
    error(' Problem in PARAFAC initialization - Not set correct')
end

% Check for signs and reflect if appropriate

for f=1:Fac
    for m=1:ord-1
        if sign(sum(Factors{m}(:,f)<0)) & FixMode(m)==0
            contin=1;
            for m2 = m+1:ord
                if contin & FixMode(m2)==0
                    if sign(sum(Factors{m2}(:,f)<0))
                        Factors{m}(:,f)=-Factors{m}(:,f);
                        Factors{m2}(:,f)=-Factors{m2}(:,f);
                        contin=0;
                    end
                end
            end
        end
    end
end
% Convert to old format


if iscell(Factors)
    ff = [];
    for f=1:length(Factors)
        ff=[ff;Factors{f}(:)];
    end
    Factors = ff;
end


% ALTERNATING LEAST SQUARES
err=SSX;
f=2*crit;
it=0;
connew=2;conold=1; % for missing values
ConstraintsNotRight = 0; % Just to ensure that iterations are not stopped if constraints are not yet fully imposed

if showfit~=-1
    disp(' ')
    if iscell(const)
        disp(' Sum-of-Squares   Iterations  Explained')
        disp(' of residuals                 variation')
    elseif any(const==4)|any(const==5)
        disp(' Sum-of-absolute   Iterations  Explained')
        disp(' residuals                     variation')
    else
        disp(' Sum-of-Squares   Iterations  Explained')
        disp(' of residuals                 variation')
    end
end


while (((f>crit) | (norm(connew-conold)/norm(conold)>MissConvCrit) | ConstraintsNotRight) & it<maxit)|~ nonneg_obeyed
    conold=connew; % for missing values
    it=it+1;
    acc=acc+1;
    
    if strcmp(accel_pattern,'paired') % this chooses Bro standard two-steps version
        if it==1
            startiter=-acc+1;
        end
        
        if acc==do_acc;
            Load_o1=Factors;
        end
        if acc==do_acc+1;
            acc=0;Load_o2=Factors;
            Factors=Load_o1+(Load_o2-Load_o1)*(it^(1/acc_pow));
            % Convert to new format
            clear ff,id1 = 0;
            for i = 1:length(DimX)
                id2 = sum(DimX(1:i).*Fac);ff{i} = reshape(Factors(id1+1:id2),DimX(i),Fac);id1 = id2;
            end
            model=nmodel(ff);
            model = reshape(model,DimX(1),prod(DimX(2:end)));
            %G(it)=err;
            
            if MissMeth
                connew=model(id);
                errX=X-model;
                if iscell(const)
                    nerr=sum(sum(errX(idmiss2).^2));
                elseif DoWeight==0&~(any(const==4)|any(const==5))
                    nerr=sum(sum(errX(idmiss2).^2));
                elseif any(const==4)|any(const==5)
                    nerr=sum(sum(abs(errX(idmiss2(:)))));
                else
                    nerr=sum(sum((Weights(idmiss2).*errX(idmiss2)).^2));
                end
            else
                if iscell(const)
                    nerr=sum(sum((X-model).^2));
                elseif DoWeight==0&~(any(const==4)|any(const==5))
                    nerr=sum(sum((X-model).^2));
                elseif any(const==4)|any(const==5)
                    nerr=sum(sum(abs(X(:)-model(:))));
                else
                    nerr=sum(sum((X.*Weights-model.*Weights).^2));
                end
            end
            if nerr>err
                acc_fail=acc_fail+1;
                Factors=Load_o2;
                if acc_fail==max_fail,
                    acc_pow=acc_pow+1+1;
                    acc_fail=0;
                    if showfit~=-1
                        disp(' Reducing acceleration');
                    end
                end
            else
                if MissMeth
                    X(id)=model(id);
                end
            end
        end
    elseif strcmp(accel_pattern,'shifted')% this chooses harshman's every-iter way
        if acc==do_acc-1;
            Load_o2=Factors;
            disp([' LS shifted accel will be done every time after iter ' int2str(it)])
        elseif acc==do_acc;
            Load_o1=Load_o2; % shift the one that was latest back to first
            acc=0;Load_o2=Factors; % latest factors become the newer set
            Factors=Load_o1+(Load_o2-Load_o1)*(it^(1/acc_pow));
            % Convert to new format
            clear ff,id1 = 0;
            for i = 1:length(DimX)
                id2 = sum(DimX(1:i).*Fac);ff{i} = reshape(Factors(id1+1:id2),DimX(i),Fac);id1 = id2;
            end
            model=nmodel(ff);
            model = reshape(model,DimX(1),prod(DimX(2:end)));
            %G(it)=err;
            
            if MissMeth
                connew=model(id);
                errX=X-model;
                if iscell(const)
                    nerr=sum(sum(errX(idmiss2).^2));
                elseif (DoWeight==0&~(any(const==4)|any(const==5)))
                    nerr=sum(sum(errX(idmiss2).^2));
                elseif any(const==4)|any(const==5)
                    nerr=sum(sum(abs(errX(idmiss2(:)))));
                else
                    nerr=sum(sum((Weights(idmiss2).*errX(idmiss2)).^2));
                end
            else
                if iscell(const)
                    nerr=sum(sum((X-model).^2));
                elseif(DoWeight==0&~(any(const==4)|any(const==5)))
                    nerr=sum(sum((X-model).^2));
                elseif any(const==4)|any(const==5)
                    nerr=sum(sum(abs(X(:)-model(:))));
                else
                    nerr=sum(sum((X.*Weights-model.*Weights).^2));
                end
            end
            if nerr>err
                acc_fail=acc_fail+1;
                Factors=Load_o2;
                if acc_fail==max_fail,
                    acc_pow=acc_pow+1+1;
                    acc_fail=0;
                    if showfit~=-1
                        disp(' Reducing acceleration');
                    end
                end
            else
                if MissMeth
                    X(id)=model(id);
                end
            end
        end
    elseif strcmp(accel_pattern,'none')
        % do nothing
    else
        error('accel_pattern should be ''shifted'', ''paired'', or ''none''.')
    end
    
    %   if acc==do_acc; %     Load_o1=Factors;    %   end,    %   if,     %   acc==do_acc+1;    %     acc=0;Load_o2=Factors;    %     Factors=Load_o1+(Load_o2-Load_o1)*(it^(1/acc_pow));    %     %,    %     %     Convert to new format,    %     clear ff,id1 = 0;    %     for i,    %     = 1:length(DimX),    %       id2 = sum(DimX(1:i).*Fac);ff{i} = reshape(Factors(id1+1:id2),DimX(i),Fac);id1 = id2;    %     end,    %     model=nmodel(ff);    %     model = reshape(model,DimX(1),prod(DimX(2:end)));    %     if MissMeth,    %       connew=model(id);    %       errX=X-model;    %       if DoWeight==0,    %,    %,    %       nerr=sum(sum(errX(idmiss2).^2));    %       else,    %         nerr=sum(sum((Weights(idmiss2).*errX(idmiss2)).^2));    %       end,    %     else,    %       if DoWeight==0,    %         nerr=sum(sum((X-model).^2));    %       else,    %    %         nerr=sum(sum((X.*Weights-model.*Weights).^2));    %       end,    %     end,    %     if nerr>err,    %       acc_fail=acc_fail+1;    %       Factors=Load_o2;    %       if acc_fail==max_fail,    %         acc_pow=acc_pow+1+1;    %         acc_fail=0;    %         if showfit~=-1,    %,     %         disp(' Reducing acceleration');    %         end,    %,    %         end,    %     else,    %       if MissMeth,    %         X(id)=model(id);    %       end,    %     end,    %   end
    if ~iscell(const)
        if DoWeight==0&~(any(const==4)|any(const==5))
            for ii=ord:-1:1
                if ii==ord;
                    i=1;
                else
                    i=ii+1;
                end
                idd=[i+1:ord 1:i-1];
                l_idx2=lidx(idd,:);
                dimx=DimX(idd);
                if ~FixMode(i)
                    L1=reshape(Factors(l_idx2(1,1):l_idx2(1,2)),dimx(1),Fac);
                    if ord>2
                        L2=reshape(Factors(l_idx2(2,1):l_idx2(2,2)),dimx(2),Fac);
                        Z=krb(L2,L1);
                    else
                        Z = L1;
                    end
                    for j=3:ord-1
                        L1=reshape(Factors(l_idx2(j,1):l_idx2(j,2)),dimx(j),Fac);
                        Z=krb(L1,Z);
                    end
                    ZtZ=Z'*Z;
                    ZtX=Z'*X';
                    OldLoad=reshape(Factors(lidx(i,1):lidx(i,2)),DimX(i),Fac);
                    L=pfls(ZtZ,ZtX,DimX(i),const(i),OldLoad,DoWeight,Weights);
                    %L=(pinv(ZtZ+eye*600)*ZtX)';
                    Factors(lidx(i,1):lidx(i,2))=L(:);
                end
                x=zeros(prod(DimX([1:ii-1 ii+1:ord])),DimX(ii));  % Rotate X so the current last mode is the first
                x(:)=X;
                X=x';
            end
        elseif DoWeight
            for ii=ord:-1:1
                if ii==ord;
                    i=1;
                else
                    i=ii+1;
                end
                idd=[i+1:ord 1:i-1];
                l_idx2=lidx(idd,:);
                dimx=DimX(idd);
                if ~FixMode(i)
                    L1=reshape(Factors(l_idx2(1,1):l_idx2(1,2)),dimx(1),Fac);
                    if ord>2
                        L2=reshape(Factors(l_idx2(2,1):l_idx2(2,2)),dimx(2),Fac);
                        Z=krb(L2,L1);
                    else
                        Z = L1;
                    end
                    for j=3:ord-1
                        L1=reshape(Factors(l_idx2(j,1):l_idx2(j,2)),dimx(j),Fac);
                        Z=krb(L1,Z);
                    end
                    OldLoad=reshape(Factors(lidx(i,1):lidx(i,2)),DimX(i),Fac);
                    L=pfls(Z,X,DimX(i),const(i),OldLoad,DoWeight,Weights);
                    Factors(lidx(i,1):lidx(i,2))=L(:);
                end
                x=zeros(prod(DimX([1:ii-1 ii+1:ord])),DimX(ii));
                x(:)=X;
                X=x';
                x(:)=Weights;
                Weights=x';
            end
        elseif (any(const==4)|any(const==5))
            for ii=ord:-1:1
                if ii==ord;
                    i=1;
                else
                    i=ii+1;
                end
                idd=[i+1:ord 1:i-1];
                l_idx2=lidx(idd,:);
                dimx=DimX(idd);
                if ~FixMode(i)
                    L1=reshape(Factors(l_idx2(1,1):l_idx2(1,2)),dimx(1),Fac);
                    if ord>2
                        L2=reshape(Factors(l_idx2(2,1):l_idx2(2,2)),dimx(2),Fac);
                        Z=krb(L2,L1);
                    else
                        Z = L1;
                    end
                    for j=3:ord-1
                        L1=reshape(Factors(l_idx2(j,1):l_idx2(j,2)),dimx(j),Fac);
                        Z=krb(L1,Z);
                    end
                    OldLoad=reshape(Factors(lidx(i,1):lidx(i,2)),DimX(i),Fac);
                    if any(const==5);
                        L=plae(X',Z,OldLoad',1)';
                    else
                        L=plae(X',Z,OldLoad',0)';
                    end
                    Factors(lidx(i,1):lidx(i,2))=L(:);
                end
                x=zeros(prod(DimX([1:ii-1 ii+1:ord])),DimX(ii));
                x(:)=X;
                X=x';
            end
        end
        
    else
        for ii=ord:-1:1
            if ii==ord;
                i=1;
            else
                i=ii+1;
            end
            idd=[i+1:ord 1:i-1];
            l_idx2=lidx(idd,:);
            dimx=DimX(idd);
            if ~FixMode(i)
                L1=reshape(Factors(l_idx2(1,1):l_idx2(1,2)),dimx(1),Fac);
                if ord>2
                    L2=reshape(Factors(l_idx2(2,1):l_idx2(2,2)),dimx(2),Fac);
                    Z=krb(L2,L1);
                else
                    Z = L1;
                end
                for j=3:ord-1
                    L1=reshape(Factors(l_idx2(j,1):l_idx2(j,2)),dimx(j),Fac);
                    Z=krb(L1,Z);
                end
                OldLoad=reshape(Factors(lidx(i,1):lidx(i,2)),DimX(i),Fac);
                L=compresnonneg(X',Z,OldLoad,const{i},it);
                Factors(lidx(i,1):lidx(i,2))=L(:);
            end
            x=zeros(prod(DimX([1:ii-1 ii+1:ord])),DimX(ii));  % Rotate X so the current last mode is the first
            x(:)=X;
            X=x';
        end
    end
    
    
    
    % POSTPROCES LOADINGS (ALL VARIANCE IN FIRST MODE)
    if ~any(FixMode)
        A=reshape(Factors(lidx(1,1):lidx(1,2)),DimX(1),Fac);
        for i=2:ord
            B=reshape(Factors(lidx(i,1):lidx(i,2)),DimX(i),Fac);
            for ff=1:Fac
                A(:,ff)=A(:,ff)*norm(B(:,ff));
                B(:,ff)=B(:,ff)/norm(B(:,ff));
            end
            Factors(lidx(i,1):lidx(i,2))=B(:);
        end
        Factors(lidx(1,1):lidx(1,2))=A(:);
    end
    % APPLY SIGN CONVENTION IF NO FIXED MODES
    %  FixMode=1
    if ~iscell(const)
        if ~any(FixMode)&~(any(const==2)|any(const==3))
            Sign = ones(1,Fac);
            for i=ord:-1:2
                A=reshape(Factors(lidx(i,1):lidx(i,2)),DimX(i),Fac);
                Sign2=ones(1,Fac);
                for ff=1:Fac
                    [out,sig]=max(abs(A(:,ff)));
                    Sign(ff) = Sign(ff)*sign(A(sig,ff));
                    Sign2(ff) = sign(A(sig,ff));
                end
                A=A*diag(Sign2);
                Factors(lidx(i,1):lidx(i,2))=A(:);
            end
            A=reshape(Factors(lidx(1,1):lidx(1,2)),DimX(1),Fac);
            A=A*diag(Sign);
            Factors(lidx(1,1):lidx(1,2))=A(:);
        end
    end
    % Check if nonneg_obeyed
    if ~iscell(const)
        for i=1:ord
            if const(i)==2|const(i)==3
                A=reshape(Factors(lidx(i,1):lidx(i,2)),DimX(i),Fac);
                if any(A(:))<0
                    nonneg_obeyed=0;
                end
            end
        end
    end
    
    % EVALUATE SOFAR
    % Convert to new format
    clear ff,id1 = 0;
    for i = 1:length(DimX)
        id2 = sum(DimX(1:i).*Fac);
        ff{i} = reshape(Factors(id1+1:id2),DimX(i),Fac);
        id1 = id2;
    end
    model=nmodel(ff);
    model = reshape(model,DimX(1),prod(DimX(2:end)));
    if MissMeth  % Missing values present
        connew=model(id);
        X(id)=model(id);
        errold=err;
        errX=X-model;
        if iscell(const)
            err=sum(sum(errX(idmiss2).^2));
        elseif (DoWeight==0&~(any(const==4)|any(const==5)))
            err=sum(sum(errX(idmiss2).^2));
        elseif any(const==4)|any(const==5)
            err=sum(sum(abs(errX(idmiss2(:)))));
        else
            err=sum(sum((Weights(idmiss2).*errX(idmiss2)).^2));
        end
    else
        errold=err;
        if iscell(const)
            err=sum(sum((X(:)-model(:)).^2));
        elseif (DoWeight==0&~(any(const==4)|any(const==5)))
            err=sum(sum((X-model).^2));
        elseif any(const==4)|any(const==5)
            err=sum(sum(abs(X(:)-model(:))));
        else
            err=sum(sum((Weights.*(X-model)).^2));
        end
    end
    if err/SSX<1000*eps, % Getting close to the machine uncertainty => stop
        disp(' WARNING')
        disp(' The misfit is approaching the machine uncertainty')
        disp(' If pure synthetic data is used this is OK, otherwise if the')
        disp(' data elements are very small it might be appropriate ')
        disp(' to multiply the whole array by a large number to increase')
        disp(' numerical stability. This will only change the solution ')
        disp(' by a scaling constant')
        f = 0;
    else
        f=abs((err-errold)/err);
        if f<crit % Convergence: then check that constraints are fulfilled
            if ~iscell(const)
                if any(const==2)|any(const==3) % If nnls or unimodality imposed
                    for i=1:ord % Extract the
                        if const(i)==2|const(i)==3 % If nnls or unimodality imposed
                            Loadd = Factors(sum(DimX(1:i-1))*Fac+1:sum(DimX(1:i))*Fac);
                            if any(Loadd<0)
                                ConstraintsNotRight=1;
                            else
                                ConstraintsNotRight=0;
                            end
                        end
                    end
                end
            end
        end
    end
    
    if it/showfit-round(it/showfit)==0
        if showfit~=-1,
            ShowPhi=ShowPhi+1;
            if ShowPhi==ShowPhiWhen,
                ShowPhi=0;
                if showfit~=-1,
                    disp(' '),
                    disp('    Tuckers congruence coefficient'),
                    % Convert to new format
                    clear ff,id1 = 0;
                    for i = 1:length(DimX)
                        id2 = sum(DimX(1:i).*Fac);ff{i} = reshape(Factors(id1+1:id2),DimX(i),Fac);id1 = id2;
                    end
                    [phi,out]=ncosine(ff,ff);
                    disp(phi),
                    if MissMeth
                        fprintf(' Change in estim. missing values %12.10f',norm(connew-conold)/norm(conold));
                        disp(' ')
                        disp(' ')
                    end
                    disp(' Sum-of-Squares   Iterations  Explained')
                    disp(' of residuals                 variation')
                end
            end
            if iscell(const)
                PercentExpl=100*(1-(sum(sum((X-model).^2))/SSX));
            elseif (DoWeight==0&~(any(const==4)|any(const==5)))
                PercentExpl=100*(1-err/SSX);
            elseif any(const==4)|any(const==5)
                PercentExpl=100*(1-sum(abs((X(:)-model(:))))/SSX);
            else
                PercentExpl=100*(1-(sum(sum((X-model).^2))/SSX));
            end
            fprintf(' %12.10f       %g        %3.4f    \n',err,it,PercentExpl);
            if Plt==2|Plt==3
                % Convert to new format
                clear ff,id1 = 0;
                for i = 1:length(DimX)
                    id2 = sum(DimX(1:i).*Fac);ff{i} = reshape(Factors(id1+1:id2),DimX(i),Fac);id1 = id2;
                end
                pfplot(reshape(X,DimX),ff,Weights',[0 0 0 0 0 0 0 1]);
                drawnow
            end
        end
    end
    
    % Make safety copy of loadings and initial parameters in temp.mat
    if it/50-round(it/50)==0
        save temp Factors
    end
    
    % JUDGE FIT
    if err>errold
        NumberOfInc=NumberOfInc+1;
    end
    % POSTPROCESS. IF PCA on two-way enforce orth in both modes.
    
end % while f>crit

if DoingPCA
    A=reshape(Factors(lidx(1,1):lidx(1,2)),DimX(1),Fac);
    B=reshape(Factors(lidx(i,1):lidx(i,2)),DimX(i),Fac);
    [u,s,v]=svd(A*B',0);
    A = u(:,1:size(A,2))*s(1:size(A,2),1:size(A,2));
    B = u(:,1:size(B,2));
    Factors = [A(:);B(:)];
end


% CALCULATE TUCKERS CONGRUENCE COEFFICIENT
if showfit~=-1 & DimX(1)>1
    disp(' '),disp('   Tuckers congruence coefficient')
    % Convert to new format
    clear ff,id1 = 0;
    for i = 1:length(DimX)
        id2 = sum(DimX(1:i).*Fac);ff{i} = reshape(Factors(id1+1:id2),DimX(i),Fac);id1 = id2;
    end
    [phi,out]=ncosine(ff,ff);
    disp(phi)
    disp(' ')
    if max(max(abs(phi)-diag(diag(phi))))>.85
        disp(' ')
        disp(' ')
        disp(' WARNING, SOME FACTORS ARE HIGHLY CORRELATED.')
        disp(' ')
        disp(' You could decrease the number of components. If this')
        disp(' does not help, try one of the following')
        disp(' ')
        disp(' - If systematic variation is still present you might')
        disp('   wanna decrease your convergence criterion and run')
        disp('   one more time using the loadings as initial guess.')
        disp(' ')
        disp(' - Or use another preprocessing (check for constant loadings)')
        disp(' ')
        disp(' - Otherwise try orthogonalising some modes,')
        disp(' ')
        disp(' - Or use Tucker3/Tucker2,')
        disp(' ')
        disp(' - Or a PARAFAC with some modes collapsed (if # modes > 3)')
        disp(' ')
    end
end


% SHOW FINAL OUTPUT
if iscell(const)
    PercentExpl=100*(1-sum(sum((X-model).^2))/SSX);
elseif DoWeight==0&~(any(const==4)|any(const==5))
    PercentExpl=100*(1-err/SSX);
elseif any(const==5)
    PercentExpl=100*(1-sum(abs(X(:)-model(:)))/SSX);
else
    PercentExpl=100*(1-sum(sum((X-model).^2))/SSX);
end
if showfit~=-1
    fprintf(' %12.10f       %g        %3.4f \n',err,it,PercentExpl);
    if NumberOfInc>0
        disp([' There were ',num2str(NumberOfInc),' iterations that increased fit']);
    end
end

% POSTPROCES LOADINGS (ALL VARIANCE IN FIRST MODE)
if Options(4)==0|Options(4)==1
    A=reshape(Factors(lidx(1,1):lidx(1,2)),DimX(1),Fac);
    for i=2:ord
        if ~FixMode(i)
            B=reshape(Factors(lidx(i,1):lidx(i,2)),DimX(i),Fac);
            for ff=1:Fac
                A(:,ff)=A(:,ff)*norm(B(:,ff));
                B(:,ff)=B(:,ff)/norm(B(:,ff));
            end
            Factors(lidx(i,1):lidx(i,2))=B(:);
        end
    end
    Factors(lidx(1,1):lidx(1,2))=A(:);
    if showfit~=-1
        disp(' ')
        disp(' Components have been normalized in all but the first mode')
    end
end

% PERMUTE SO COMPONENTS ARE IN ORDER AFTER VARIANCE DESCRIBED (AS IN PCA) IF NO FIXED MODES
if ~any(FixMode)
    A=reshape(Factors(lidx(1,1):lidx(1,2)),DimX(1),Fac);
    [out,order]=sort(diag(A'*A));
    order=flipud(order);
    A=A(:,order);
    Factors(lidx(1,1):lidx(1,2))=A(:);
    for i=2:ord
        B=reshape(Factors(lidx(i,1):lidx(i,2)),DimX(i),Fac);
        B=B(:,order);
        Factors(lidx(i,1):lidx(i,2))=B(:);
    end
    if showfit~=-1
        disp(' Components have been ordered according to contribution')
    end
elseif showfit ~= -1
    disp(' Some modes fixed hence no sorting of components performed')
end


% TOOLS FOR JUDGING SOLUTION
if nargout>3
    x=X;
    if MissMeth
        x(id)=NaN*id;
    end
    % Convert to new format
    clear ff,id1 = 0;
    for i = 1:length(DimX)
        id2 = sum(DimX(1:i).*Fac);ff{i} = reshape(Factors(id1+1:id2),DimX(i),Fac);id1 = id2;
    end
    any(isnan(x(:)))
    any(isnan(ff{1}(:)))
    any(isnan(ff{2}(:)))
    any(isnan(ff{3}(:)))
    corcondia=corcond(reshape(x,DimX),ff,Weights,0);
end

% APPLY SIGN CONVENTION IF NO FIXED MODES
%  FixMode=1
if ~iscell(const)
    if ~any(FixMode)&~(any(const==2)|any(const==3))
        Sign = ones(1,Fac);
        for i=ord:-1:2
            A=reshape(Factors(lidx(i,1):lidx(i,2)),DimX(i),Fac);
            Sign2=ones(1,Fac);
            for ff=1:Fac
                [out,sig]=max(abs(A(:,ff)));
                Sign(ff) = Sign(ff)*sign(A(sig,ff));
                Sign2(ff) = sign(A(sig,ff));
            end
            A=A*diag(Sign2);
            Factors(lidx(i,1):lidx(i,2))=A(:);
        end
        A=reshape(Factors(lidx(1,1):lidx(1,2)),DimX(1),Fac);
        A=A*diag(Sign);
        Factors(lidx(1,1):lidx(1,2))=A(:);
        
        %   % Instead of above, do signs so as to make them as "natural" as possible
        %   Factors = signswtch(Factors,reshape(X,DimX));
        %   DIDN't WORK (TOOK AGES FOR 7WAY DATA)
        
        
        if showfit~=-1
            disp(' Components have been reflected according to convention')
        end
    end
end

% Convert to new format
clear ff,id1 = 0;
for i = 1:length(DimX)
    id2 = sum(DimX(1:i).*Fac);ff{i} = reshape(Factors(id1+1:id2),DimX(i),Fac);id1 = id2;
end
Factors = ff;


if Plt==1|Plt==2|Plt==3
    %   if Fac<6&Plt~=3&order>2&ord>2
    if Fac<6&Plt~=3&ord>2
        pfplot(reshape(X,DimX),ff,Weights,ones(1,8));
    else
        pfplot(reshape(X,DimX),ff,Weights,[1 1 0 1 1 1 1 1]);
        if ord>2
            disp(' Core consistency plot not shown because it requires large memory')
            disp(' It can be made writing pfplot(X,Factors,[Weights],[0 0 1 0 0 0 0 0]');
        else
            disp(' Core consistency not applicable for two-way data')
        end
    end
end

% Show which criterion stopped the algorithm
if showfit~=-1
    if ((f<crit) & (norm(connew-conold)/norm(conold)<MissConvCrit))
        disp(' The algorithm converged')
    elseif it==maxit
        disp(' The algorithm did not converge but stopped because the')
        disp(' maximum number of iterations was reached')
    elseif f<eps
        disp(' The algorithm stopped because the change in fit is now')
        disp(' smaller than the machine uncertainty.')
    else
        disp(' Algorithm stopped for some mysterious reason')
    end
end

function swloads = signswtch(loads,X);

%SIGNSWTCH switches sign of multilinear models so that signs are in
%accordance with majority of data
%
%
% I/O swloads = signswtch(loads,X);
%
% Factors must be a cell with the loadings. If Tucker or NPLS, then the
% last element of the cell must be the core array

try % Does not work in older versions of matlab
    warning('off','MATLAB:divideByZero');
end
sizeX=size(X);
order = length(sizeX);
for i=1:order;
    F(i) = size(loads{i},2);
end


if isa(X,'dataset')% Then it's a SDO
    inc=X.includ;
    X = X.data(inc{:});
end


% Compare centered X with center loading vector

if length(loads)==order % PARAFAC
    % go through each component and then update in the end
    for m = 1:order % For each mode determine the right sign
        for f=1:F(1) % one factor at the time
            s=[];
            a = loads{m}(:,f);
            x = permute(X,[m 1:m-1 m+1:order]);
            for i=1:size(x(:,:),2); % For each column
                id = find(~isnan(x(:,i)));
                if length(id)>1
                    try
                        c = corrcoef(x(id,i),a(id));
                    catch
                        disp('Oops - something wrong in signswtch - please send a note to rb@kvl.dk')
                        whos
                    end
                    if isnan(c(2,1))
                        s(i)=0;
                    else
                        s(i) = c(2,1)*length(id); % Weigh correlation by number of elements so many-miss columns don't influence too much
                    end
                else
                    s(i) = 0;
                end
            end
            S(m,f) = sum(s);
        end
    end
    
    % Use S to switch signs. If the signs of S (for each f) multiply to a
    % positive number the switches are performed. If not, the mode of the
    % negative one with the smallest absolute value is not switched.
    
    for f = 1:F(1)
        if sign(prod(S(:,f)))<1 % Problem: make the smallest negative positive to avoid switch of that
            id = find(S(:,f)<0);
            [a,b]=min(abs(S(id,f)));
            S(id(b(1)),f)=-S(id(b(1)),f);
        end
    end
    % Now ok, so switch what needs to be switched
    for f = 1:F(1)
        for m = 1:order
            if sign(S(m,f))<1
                loads{m}(:,f)=-loads{m}(:,f);
            end
        end
    end
    
elseif length(loads)==(order+1) % NPLS/Tucker
    
    % go through each mode and update and correct core accordinglu
    for m = 1:order % For each mode determine the right sign
        for f=1:F(m) % one factor at the time
            a = loads{m}(:,f);
            x = permute(X,[m 1:m-1 m+1:order]);
            for i=1:size(x(:,:),2); % For each column
                id = find(~isnan(x(:,i)));
                if length(id)>1
                    c = corrcoef(x(id,i),a(id));
                    if isnan(c(2,1))
                        s(i)=0;
                    else
                        s(i) = c(2,1)*length(id); % Weigh correlation by number of elements so many-miss columns don't influence too much
                    end
                else
                    s(i) = 0;
                end
            end
            if sum(s) < 0
                % turn around
                loads{m}(:,f) = -loads{m}(:,f);
                
                % Then switch the core accordingly
                G = loads{order+1};
                G = permute(G,[m 1:m-1 m+1:order]);
                sizeG = size(G);
                G = reshape(G,sizeG(1),prod(sizeG)/sizeG(1));
                G(f,:) = -G(f,:);
                G = reshape(G,sizeG);
                G = ipermute(G,[m 1:m-1 m+1:order]);
                loads{order+1} = G;
            end
        end
    end
    
    
else
    error('Unknown model type in SIGNS.M')
end

swloads = loads;

function S = plae(X,A,Sold,nonnegativity);

% Least Absolute Error ``pseudoinverse'':
% iteratively solves the problem of min wrt S sum(sum(abs(X-A*S)))
% monotonically convergent in terms of LAE,
% but not guaranteed in general to find globally opt soln.
% ALL parameters are real-valued
% uses wmf.m
% complexity per complete iteration is order of F*N*I*logI
% where S is FxN, and A is IxF
% Assumes A is tall and full rank
% N. Sidiropoulos, April 12, 2000
% RB 2005, Added nonnegativity

[I,F]=size(A);
[I,N]=size(X);

if nargin<3
    S = pinv(A)*X; % use LS-solution as good initialization
else
    S = Sold;
end
LAE = sum(sum(abs(X-A*S)));
SMALLNUMBER = eps*10;
MAXNUMITER = 30;

%fprintf('LAE = %12.10f\n',LAE);
LAEold = 2*LAE;
LAEinit = LAE;
it     = 0;

while abs((LAE-LAEold)/LAEold) > SMALLNUMBER & it < MAXNUMITER & LAE > 10*eps
    it=it+1;
    LAEold=LAE;
    
    % update elements of S one by one:
    
    for f=1:F,
        Y = X - (A*S - A(:,f)*S(f,:));
        for n=1:N,
            S(f,n) = wmf(Y(:,n)./A(:,f),A(:,f));
        end
        if nonnegativity
            S(f,find(S(f,:)<0)) = 0;
        end
    end
    
    % compute new LAE:
    
    LAE = sum(sum(abs(X-A*S)));
    %fprintf('LAE = %12.10f\n',LAE);
    
end % while loop

function a = wmf(x,w);

% weighted median filter with non-negative real weights
% input (vector x) and output (scalar a) are both real-valued
% minimize sum_i abs(w(i)) * abs(x(i) - a)
% N. Sidiropoulos, April 12, 2000
% Ref: cf. e.g., Yang et al, IEEE Trans. Signal Proc. 43(3):591-592, Mar. 1995
% Complexity is NlogN, N=length(x), due to the sorting operation

[s,p] = sort(x);
absw = abs(w);
sw = absw(p);
t = 0.5*sum(absw);
N=length(x);
psum=0;
for n=N:-1:1,
    psum = psum + sw(n);
    if (psum >= t)
        a = s(n);
        break;
    end
end


function B=compresnonneg(X,Z,Bold,V,it)

% Find B minimizing ||X-Z*B'|| subject to V*B>=0
% B=compresnonneg(X,Z,V);

[I,J]=size(X);
[Iz,F]=size(Z);
[BJ,~]=size(Bold);
[VI,VJ]=size(V);
z=sparse(I*J,F*J);
for i=1:J;
    z((i-1)*Iz+1:i*Iz,(i-1)*F+1:i*F)=Z;
end

VV=sparse(VI*BJ,VJ*F);
co=0;
for i=1:VI;
    for f=1:F
        co=co+1;
        VV(co,f:F:end)=V(i,:);
    end
end

if it<10 % Impose constraints slowly during first ten iteration
    mx = -max(X(:));
    mx=mx*(10/(10-it));
    thrs=repmat(mx,size(VV,1),1);
else
    thrs=sparse(size(VV,1),1);
end
B = lsi2(z'*z,z'*X(:),VV,thrs);
B=reshape(B,F,BJ)';

function x = lsi2(EE,Ef,G,h);
% x = lsi2(EE,Ef,G,h)
% Least squares with inequality constraints.
% Solves the following problem:
%       min || Ex - f || subject to Gx >= h
% Here input: EE = E'E and Ef=E'f

R=chol(EE);
RR=inv(R);
EEinv=RR*RR';
GG=G*RR;
hh=h-G*EEinv*Ef;
z=ldp(GG,hh);
if ~any(isnan(z))
    x=RR*z+EEinv*Ef;
else
    x=NaN;
end

function b=ldp(G,h);

% find min||b|| subject to Gb=>h
%
% From Lawson & Hanson 74, p. 165
%
% Copyright 1998
% Rasmus Bro
% rasmus@optimax.dk/rb@kvl.dk


[I,J]=size(G);
E=[G';h'];
f=[zeros(J,1);1];
uhat=fastnnls(E'*E,E'*f);
r=E*uhat-f;

if norm(r)==0
    disp(' No solution to LDP problem')
elseif r(J+1)~=0
    b=-r(1:J)/r(J+1);
else
    b=NaN;
end


function [Consistency,G,stdG,Target]=corcond(X,Factors,Weights,Plot);

%CORCOND Core consistency for PARAFAC model
%
% See also:
% 'unimodal' 'monreg' 'fastnnls'
%
% CORe CONsistency DIAgnostics (corcondia)
% Performs corcondia of a PARAFAC model and returns the cocote plot
% as well as the degree of consistency (100 % is max).
%
% Consistency=corcond(X,Factors,Weights,Plot);
%
% INPUT
% X        : Data array
% Factors  : Factors given in standard format as a cell array
% Weights  : Optional weights (otherwise skip input or give an empty array [])
% Plot     = 0 or not given => no plots are produced
%          = 1              => normal corcondia plot
%          = 2              => corcondia plot with standard deviations
%
% OUTPUT
% The core consistency given as the percentage of variation in a Tucker3 core
% array consistent with the theoretical superidentity array. Max value is 100%
% Consistencies well below 70-90% indicates that either too many components
% are used or the model is otherwise mis-specified.
%

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.

% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $
% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $
% $ Version 2.01 $ Feb 2003 $ replaced regg with t3core when weights are used $ RB $ Not compiled $

DimX = size(X);
X = reshape(X,DimX(1),prod(DimX(2:end)));
Fac = size(Factors{1},2);

if nargin<4
    Plot=0;
end
if nargin<3
    Weights=0;
end

ord=length(DimX);
l_idx=0;
for i=1:ord
    l_idx=[l_idx sum(DimX(1:i))*Fac];
end


% Scale all loadings to same magnitude
magn=ones(Fac,1);
for i=1:ord
    L=Factors{i};
    for f=1:Fac
        magn(f)=magn(f)*norm(L(:,f));
        L(:,f)=L(:,f)/norm(L(:,f));
    end
    Factors{i}=L;
end
% Magn holds the singular value of each component. Scale each loading vector by
% the cubic root (if three-way) so all loadings of a component have the same variance

magn = magn.^(1/ord);
for i=1:ord
    L=Factors{i};
    for f=1:Fac
        L(:,f)=L(:,f)*magn(f);
    end
    Factors{i}=L;
end


% Make diagonal array holding the magnitudes
Ident=nident(Fac,ord);
if Fac>1
    DimIdent=ones(1,ord)*Fac;
    Ident=nshape(reshape(Ident,DimIdent),ord);
end

% Make matrix of Kronecker product of all loadings expect the large; Z = kron(C,B ... )
NewFac=[];
NewFacNo=[];
for i=ord:-1:1
    Z=Factors{i};
    % Check its of full rank or adjust core and use less columns
    rankZ=rank(Z);
    if rankZ<Fac
        %OLD out=Z(:,rankZ+1:Fac);Z=Z(:,1:rankZ);H=[[eye(rankZ)] pinv(Z)*out];Ident=H*Ident;
        [q,r]=qr(Z);
        Ident=r*Ident;
        Z=q;
        DimIdent(i)=size(r,1);
    end
    if i>1&Fac>1
        Ident=nshape(reshape(Ident,DimIdent([i:ord 1:i-1])),ord);
    end
    NewFac{i}=Z;
    NewFacNo=[rankZ NewFacNo];
end
Factors=NewFac;
Fac=NewFacNo;
if nargin<3
    [G,stdG]=regg(reshape(X,DimX),Factors,Weights); %Doesn't work with weights
else
    G=t3core(reshape(X,DimX),Factors,Weights);
    stdG = G; % Arbitrary (not used)
end

DimG = size(G);
G = G(:);

Ident=Ident(:);
Target=Ident;
[a,b]=sort(abs(Ident));
b=flipud(b);
Ident=Ident(b);
GG=G(b);
stdGG=stdG(b);
bNonZero=find(Ident);
bZero=find(~Ident);

ssG=sum(G(:).^2);
Consistency=100*(1-sum((Target-G).^2)/ssG);


if Plot
    clf
    Ver=version;
    Ver=Ver(1);
    if Fac>1
        eval(['set(gcf,''Name'',''Diagonality test'');']);
        if Ver>4
            plot([Ident(bNonZero);Ident(bZero)],'y','LineWidth',3)
            hold on
            plot(GG(bNonZero),'ro','LineWidth',3)
            plot(length(bNonZero)+1:prod(Fac),GG(bZero),'gx','LineWidth',3)
            if Plot==2
                line([[1:length(G)];[1:length(G)]],[GG GG+stdGG]','LineWidth',1,'Color',[0 0 0])
                line([[1:length(G)];[1:length(G)]],[GG GG-stdGG]','LineWidth',1,'Color',[0 0 0])
            end
            hold off
            title(['Core consistency ',num2str(Consistency),'% (yellow target)'],'FontWeight','bold','FontSize',12)
        else
            plot([Ident(bNonZero);Ident(bZero)],'y')
            hold on
            plot(GG(bNonZero),'ro')
            plot(length(bNonZero)+1:prod(Fac),GG(bZero),'gx')
            if Plot==2
                line([[1:length(G)];[1:length(G)]],[GG GG+stdGG]','LineWidth',1,'Color',[0 0 1])
                line([[1:length(G)];[1:length(G)]],[GG GG-stdGG]','LineWidth',1,'Color',[0 0 1])
            end
            hold off
            title(['Core consistency ',num2str(Consistency),'% (yellow target)'])
        end
        xlabel('Core elements (green should be zero/red non-zero)')
        ylabel('Core Size')
    else
        eval(['set(gcf,''Name'',''Diagonality test'');']);
        title(['Core consistency ',num2str(Consistency),'% (yellow target)'])
        xlabel('Core elements (green should be zero/red non-zero)')
        ylabel('Size')
        plot(GG(bNonZero),'ro')
        title(['Core consistency ',num2str(Consistency),'%'])
        xlabel('Core elements (red non-zero)')
        ylabel('Core Size')
    end
end

G = reshape(G,DimG);

function [G,stdG]=regg(X,Factors,Weights);

%REGG Calculate Tucker core
%
% Calculate Tucker3 core

% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $

DimX = size(X);
X = reshape(X,DimX(1),prod(DimX(2:end)));
Fac = size(Factors{1},2);

ord=length(DimX);
if ord<3
    disp(' ')
    disp(' !!Corcondia only applicable for three- and higher-way arrays!!')
    return
end

if length(Fac)==1
    for i=1:length(Factors)
        Fac(i) = size(Factors{i},2);
    end
end
vecX=X(:); % Vectorize X

% Make sure Weights are defined (as ones if none given)
if nargin<3
    Weights=ones(size(X));
end
if length(Weights(:))~=length(X(:));
    Weights=ones(size(X));
end
Weights=Weights(:);

% Set weights of missing elements to zero
id=find(isnan(vecX));
Weights(id)=zeros(size(id));
vecX(id)=zeros(size(id));

% Create Kronecker product of all but the last mode loadings
L2 = Factors{end-1};
L1 = Factors{end-2};
Z = kron(L2,L1);
for o=ord-3:-1:1
    Z = kron(Z,Factors{o});
end


% Make last mode loadings, L
L=Factors{end};

% We want to fit the model ||vecX - Y*vecG||, where Y = kron(L,Z), but
% we calculate Y'Y and Y'vecX by summing over k
J=prod(DimX(1:ord-1));
Ytx = 0;
YtY = 0;
for k=1:DimX(ord)
    W=Weights((k-1)*J+1:k*J);
    WW=(W.^2*ones(1,prod(Fac)));
    Yk  = kron(L(k,:),Z);
    Ytx = Ytx + Yk'*(W.*vecX((k-1)*J+1:k*J));
    YtY = YtY + (Yk.*WW)'*Yk;
end

G=pinv(YtY)*Ytx;

if nargout>1
    se = (sum(vecX.^2) + G'*YtY*G -G'*Ytx);
    mse = se/(length(vecX)-length(G));
    stdG=sqrt(diag(pinv(YtY))*mse);
end
G = reshape(G,Fac);

function [A,B,C,fit]=dtld(X,F,SmallMode);

%DTLD direct trilinear decomposition
%
% See also:
% 'gram', 'parafac'
%
%
% DIRECT TRILINEAR DECOMPOSITION
%
% calculate the parameters of the three-
% way PARAFAC model directly. The model
% is not the least-squares but will be close
% to for precise data with little model-error
%
% This implementation works with an optimal
% compression using least-squares Tucker3 fitting
% to generate two pseudo-observation matrices that
% maximally span the variation of all samples. per
% default the mode of smallest dimension is compressed
% to two samples, while the remaining modes are
% compressed to dimension F.
%
% For large arrays it is fastest to have the smallest
% dimension in the first mode
%
% INPUT
% [A,B,C]=dtld(X,F);
% X is the I x J x K array
% F is the number of factors to fit
% An optional parameter may be given to enforce which
% mode is to be compressed to dimension two
%

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.


% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $
% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $
% $ Version 1.03 $ Date 25. April 1999 $ Not compiled $

DimX = size(X);
X = reshape(X,DimX(1),prod(DimX(2:end)));

DontShowOutput = 1;

%rearrange X so smallest dimension is in first mode


if nargin<4
    [a,SmallMode] = min(DimX);
    X = nshape(reshape(X,DimX),SmallMode);
    DimX = DimX([SmallMode 1:SmallMode-1 SmallMode+1:3]);
    Fac   = [2 F F];
else
    X = nshape(reshape(X,DimX),SmallMode);
    DimX = DimX([SmallMode 1:SmallMode-1 SmallMode+1:3]);
    Fac   = [2 F F];
end
f=F;
if F==1;
    Fac   = [2 2 2];
    f=2;
end


if DimX(1) < 2
    error(' The smallest dimension must be > 1')
end

if any(DimX(2:3)-Fac(2:3)<0)
    error(' This algorithm requires that two modes are of dimension not less the number of components')
end



% Compress data into a 2 x F x F array. Only 10 iterations are used since exact SL fit is insignificant; only obtaining good truncated bases is important
[Factors,Gt]=tucker(reshape(X,DimX),Fac,[0 0 0 0 NaN 10]);
% Convert to old format
Gt = reshape(Gt,size(Gt,1),prod(size(Gt))/size(Gt,1));

[At,Bt,Ct]=fac2let(Factors);

% Fit GRAM to compressed data
[Bg,Cg,Ag]=gram(reshape(Gt(1,:),f,f),reshape(Gt(2,:),f,f),F);

% De-compress data and find A


BB = Bt*Bg;
CC = Ct*Cg;
AA = X*pinv(krb(CC,BB)).';

if SmallMode == 1
    A=AA;
    B=BB;
    C=CC;
elseif SmallMode == 2
    A=BB;
    B=AA;
    C=CC;
elseif SmallMode == 3
    A=BB;
    B=CC;
    C=AA;
end

fit = sum(sum(abs(X - AA*krb(CC,BB).').^2));
if ~DontShowOutput
    disp([' DTLD fitted raw data with a sum-squared error of ',num2str(fit)])
end

function [Factors]=ini(X,Fac,MthFl,IgnFl)
%INI initialization of loadings
%
% function [Factors]=ini(X,Fac,MthFl,IgnFl)
%
% This algorithm requires access to:
% 'gsm' 'fnipals' 'missmult'
%
% ---------------------------------------------------------
%                    Initialize Factors
% ---------------------------------------------------------
%
% [Factors]=ini(X,Fac,MthFl,IgnFl);
% [Factors]=ini(X,Fac,MthFl);
%
% X        : The multi-way data.
% Fac      : Vector describing the number of factors
%            in each of the N modes.
% MthFl    : Method flag indicating what kind of
%            factors you want to initiate Factors with:
%            '1' : Random values, orthogonal
%            '2' : Normalized singular vectors, orthogonal
% IgnFl    : This feature is only valid with MthFl==2.
%            If specified, these mode(s) will be ignored,
%            e.g. IgnFl=[1 5] or IgnFl=[3] will
%            respectively not initialize modes one and
%            five, and mode three.
% Factors  : Contains, no matter what method, orthonormal
%            factors. This is the best general approach to
%            avoid correlated, hence ill-posed, problems.
%
% Note that it IS possible to initialize the factors to have
% more columns than rows, since this may be required by some
% PARAFAC models. If this is required, the 'superfluos'
% columns will be random and orthogonal columns.
% This algorithm automatically arranges the sequence of the
% initialization to minimize time and memory consumption.
% Note, if you get a warning from NIPALS about convergence has
% not been reached, you can simply ignore this. With regards
% to initialization this is not important as long as the
% factors being returned are in the range of the eigensolutions.

% $ Version 1.02 $ Date 30 Aug 1999 $ Not compiled $
% $ Version 1.0201 $ Date 21 Jan 2000 $ Not compiled $ RB removed orth of additional columns
% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.

format long
format compact

DimX = size(X);
X = reshape(X,DimX(1),prod(DimX(2:end)));

% Assign intermediaries
Show=0;
rand('seed',sum(100*clock));
MissingExist=any(isnan(X(:)));

% Initialize system variables
N=size(Fac,2);
if N==1,
    Fac=Fac*ones(1,size(DimX,2));
end;
N=size(Fac,2);

FIdx0=zeros(1,N);
FIdx1=zeros(1,N);
latest=1;
for c=1:N,
    if Fac(c)==-1,
        FIdx0(c)=0;
    else
        FIdx0(c)=latest;
        latest=latest+Fac(c)*DimX(c);
        FIdx1(c)=latest-1;
    end;
end;

% Check inputs
if ~exist('IgnFl'),
    IgnFl=[0];
end;

%Random values
if MthFl==1,
    for c=1:N,
        A=orth(rand( DimX(c) , min([Fac(c) DimX(c)]) ));
        %B=[A orth(rand(DimX(c),Fac(c)-DimX(c)))];
        B=[A rand(DimX(c),Fac(c)-DimX(c))];
        Factors(FIdx0(c):FIdx1(c))=B(:)';
    end;
    if Show>=1,
        fprintf('ini.m : Initialized using random values.\n');
    end;
else
    %Singular vectors
    Factors=rand(1,sum(~(Fac==-1).*DimX.*Fac));
    if MthFl==2 | MthFl==3
        [A Order]=sort(Fac);
        RedData=X;
        CurDimX=DimX;
        for k=1:N,
            c=Order(k);
            if Fac(c)>0,
                for c1=1:c-1;
                    newi=CurDimX(c1+1);
                    newj=prod(CurDimX)/CurDimX(c1+1);
                    RedData=reshape(RedData',newi,newj);
                end;
                Op=0;
                if MissingExist | (Op==0 & Fac(c)<=5 & (50<min(size(RedData)) & min(size(RedData))<=120)),
                    %Need to apply NIPALS
                    t0=clock;
                    A=reshape(Factors(FIdx0(c):FIdx1(c)),DimX(c),Fac(c));
                    if MissingExist
                        MissIdx=find(isnan(RedData));
                        [A,P]=fnipals(RedData,min([Fac(c) DimX(c)]),A);
                        Xm=A*P';
                        RedData(MissIdx)=Xm(MissIdx);
                        MissingExist=0;
                    else
                        [A]=fnipals(RedData,min([Fac(c) DimX(c)]),A);
                    end;
                    B=[A orth(rand(DimX(c),Fac(c)-DimX(c)))];
                    Factors(FIdx0(c):FIdx1(c))=B(:)';
                    t1=clock;
                    if Show>=2,
                        disp(['ini.m: NIPALS used ' num2str(etime(t1,t0)) ' secs. on mode ' int2str(c)]),
                    end;
                    Op=1;
                end;
                if Op==0 & (120<min(size(RedData)) & min(size(RedData))<Inf),
                    %Need to apply Gram-Schmidt
                    t0=clock;
                    C=RedData*RedData';
                    A=reshape(Factors(FIdx0(c):FIdx1(c)),DimX(c),Fac(c));
                    for i=1:3,
                        A=gsm(C*A);
                    end;
                    B=[A orth(rand(DimX(c),Fac(c)-DimX(c)))];
                    Factors(FIdx0(c):FIdx1(c))=B(:)';
                    t1=clock;
                    if Show>=2,
                        disp(['ini.m: GS used ' num2str(etime(t1,t0)) ' secs. on mode ' int2str(c)]),
                    end;
                    Op=1;
                end;
                if Op==0 & (0<min(size(RedData)) & min(size(RedData))<=200),
                    %Small enough to apply SVD
                    t0=clock;
                    if max(size(RedData))<1000
                        [U S A]=svd(RedData',0);
                    else
                        [U S A]=svds(RedData');
                    end
                    A=A(:,1:min(size(A,2),min([Fac(c) DimX(c)])));
                    if size(A,2)<Fac(c)
                        A = [A rand(size(A,1),Fac(c)-size(A,2))];
                    end
                    n_ = Fac(c)- min([Fac(c) DimX(c)]);
                    if n_>0,
                        a = rand(DimX(c),n_);
                        if DimX(c)>=n_
                            a = orth(a);
                        else
                            a = orth(a')';
                        end;
                        B=[A a];
                    else
                        Factors(FIdx0(c):FIdx1(c))=A(:)';
                    end;
                    
                    t1=clock;
                    if Show>=2,
                        disp(['ini.m: SVD used ' num2str(etime(t1,t0)) ' secs. on mode ' int2str(c)]),
                    end;
                    Op=1;
                end;
                CurDimX(c)=min([Fac(c) DimX(c)]);
                if MissingExist,
                    RedData=missmult(A',RedData);
                else
                    RedData=A'*RedData;
                end;
                %Examine if re-ordering is necessary
                if c~=1,
                    for c1=c:N,
                        if c1~=N,
                            newi=CurDimX(c1+1);
                            newj=prod(CurDimX)/newi;
                        else
                            newi=CurDimX(1);
                            newj=prod(CurDimX)/newi;
                        end;
                        RedData=reshape(RedData',newi,newj);
                    end;
                end;
            end;
        end;
        if Show>=1,
            fprintf('ini.m : Initialized using SVD and projection.\n');
        end;
    end;
end,
format
% Convert to new format
clear ff,id1 = 0;
for i = 1:length(DimX)
    id2 = sum(DimX(1:i).*Fac(1:i));
    ff{i} = reshape(Factors(id1+1:id2),DimX(i),Fac(i));id1 = id2;
end
Factors = ff;


function AB = krb(A,B);
%KRB Khatri-Rao product
%
% The columnwise Khatri-Rao-Bro product (Harshman, J.Chemom., 2002, 198-205)
% For two matrices with similar column dimension the khatri-Rao-Bro product
% is krb(A,B) = [kron(A(:,1),B(:,1)) .... kron(A(:,F),B(:,F))]
%
% I/O AB = krb(A,B);
%

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.

% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $
% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $
% $ Version 2.01 $ May 2001 $ Error in helpfile - A and B reversed $ RB $ Not compiled $

[I,F]=size(A);
[J,F1]=size(B);

if F~=F1
    error(' Error in krb.m - The matrices must have the same number of columns')
end

AB=zeros(I*J,F);
for f=1:F
    ab=B(:,f)*A(:,f).';
    AB(:,f)=ab(:);
end

function [MultPhi,Phis] = ncosine(factor1,factor2);

%NCOSINE multiple cosine/Tuckers congruence coefficient
%
% [MultPhi,Phis] = ncosine(factor1,factor2,DimX,Fac);
%
% ----------------------INPUT---------------------
%
% factor1   = cell array with loadings of one model
% factor2   = cell array with loadings of one (other) model
%     If factor1 and factor2 are identical then
%        the multiple cosine of a given solution is
%          estimated; otherwise the similarity of the
%          two different solutions is given
%
% ----------------------OUTPUT---------------------
%
% MultPhi   Is the multiple cosine of the model
% Phis      Is the cosine between components in
%          individual component matrices arranged
%          as [PhiA;PhiB ...]

% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $
% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $
% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.


% Convert to old format
Fac = size(factor1,2);
for i = 1:length(factor1)
    DimX(i) = size(factor1{i},1);
end

ff = [];
for f=1:length(factor1)
    ff=[ff;factor1{f}(:)];
end
factor1 = ff;

ff = [];
for f=1:length(factor2)
    ff=[ff;factor2{f}(:)];
end
factor2 = ff;


if length(factor1)~=length(factor2)
    error(' factor1 and factor2 must hold components of same sizes in NCOSINE.M')
end
ord=length(DimX);
l_idx=0;
Fac=length(factor1)/sum(DimX);
for o=1:ord
    l_idx=[l_idx sum(DimX(1:o))*Fac];
end
L1=reshape(factor1(1:DimX(1)*Fac),DimX(1),Fac);
L2=reshape(factor2(1:DimX(1)*Fac),DimX(1),Fac);
for f=1:Fac
    L1(:,f)=L1(:,f)/norm(L1(:,f));
    L2(:,f)=L2(:,f)/norm(L2(:,f));
end
%GT correction
Phis=L1'*L2;
%Previously: Phis=L2'*L2;
%End GT correction
MultPhi=Phis;

for i=2:ord
    L1=reshape(factor1(l_idx(i)+1:l_idx(i+1)),DimX(i),Fac);
    L2=reshape(factor2(l_idx(i)+1:l_idx(i+1)),DimX(i),Fac);
    for f=1:Fac
        L1(:,f)=L1(:,f)/norm(L1(:,f));
        L2(:,f)=L2(:,f)/norm(L2(:,f));
    end
    phi=(L1'*L2);
    MultPhi=MultPhi.*phi;
    Phis=[Phis;phi];
end



function [Xm]=nmodel(Factors,G,Om);

%NMODEL make model of data from loadings
%
% function [Xm]=nmodel(Factors,G,Om);
%
% This algorithm requires access to:
% 'neye.m'
%
%
% [Xm]=nmodel(Factors,G,Om);
%
% Factors  : The factors in a cell array. Use any factors from
%            any model.
% G        : The core array. If 'G' is not defined it is assumed
%            that a PARAFAC model is being established.
%            Use G = [] in the PARAFAC case.
% Om       : Oblique mode.
%            'Om'=[] or 'Om'=0, means that orthogonal
%                   projections are requsted. (default)
%            'Om'=1 means that the factors are oblique.
%            'Om'=2 means that the ortho/oblique is solved automatically.
%                   This takes a little additional time.
% Xm       : The model of X.
%
% Using the factors as they are (and the core, if defined) the general N-way model
% is calculated.


% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $
% $ Version 1.02 $ Date 17. Apr 1999 $ Not compiled $
% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.


for i = 1:length(Factors);
    DimX(i)=size(Factors{i},1);
end
i = find(DimX==0);
for j = 1:length(i)
    DimX(i(j)) = size(G,i(j));
end



if nargin<2, %Must be PARAFAC
    Fac=size(Factors{1},2);
    G=[];
else
    for f = 1:length(Factors)
        if isempty(Factors{f})
            Fac(f) = -1;
        else
            Fac(f) = size(Factors{f},2);
        end;
    end
end

if ~exist('Om')
    Om=[];
end;

if isempty(Om)
    Om=0;
end;

if size(Fac,2)==1,
    Fac=Fac(1)*ones(1,size(DimX,2));
end;
N=size(Fac,2);

if size(DimX,2)>size(Fac,2),
    Fac=Fac*ones(1,size(DimX,2));
end;
N=size(Fac,2);

Fac_orig=Fac;
i=find(Fac==-1);
if ~isempty(i)
    Fac(i)=zeros(1,length(i));
    Fac_ones(i)=ones(1,length(i));
end;
DimG=Fac;
i=find(DimG==0);
DimG(i)=DimX(i);

if isempty(G),
    G=neye(DimG);
end;
G = reshape(G,size(G,1),prod(size(G))/size(G,1));

% reshape factors to old format
ff = [];
for f=1:length(Factors)
    ff=[ff;Factors{f}(:)];
end
Factors = ff;


if DimG(1)~=size(G,1) | prod(DimG(2:N))~=size(G,2),
    
    help nmodel
    
    fprintf('nmodel.m   : ERROR IN INPUT ARGUMENTS.\n');
    fprintf('             Dimension mismatch between ''Fac'' and ''G''.\n\n');
    fprintf('Check this : The dimensions of ''G'' must correspond to the dimensions of ''Fac''.\n');
    fprintf('             If a PARAFAC model is established, use ''[]'' for G.\n\n');
    fprintf('             Try to reproduce the error and request help at rb@kvl.dk\n');
    return;
end;

if sum(DimX.*Fac) ~= length(Factors),
    help nmodel
    fprintf('nmodel.m   : ERROR IN INPUT ARGUMENTS.\n');
    fprintf('             Dimension mismatch between the number of elements in ''Factors'' and ''DimX'' and ''Fac''.\n\n');
    fprintf('Check this : The dimensions of ''Factors'' must correspond to the dimensions of ''DimX'' and ''Fac''.\n');
    fprintf('             You may be using results from different models, or\n');
    fprintf('             You may have changed one or more elements in ''Fac'' or ''DimX'' after ''Factors'' have been calculated.\n\n');
    fprintf('             Read the information above for information on arguments.\n');
    return;
end;

FIdx0=cumsum([1 DimX(1:N-1).*Fac(1:N-1)]);
FIdx1=cumsum([DimX.*Fac]);

if Om==0,
    Orthomode=1;
end;

if Om==1,
    Orthomode=0;
end;

if Om==2,
    Orthomode=1;
    for c=1:N,
        if Fac_orig(c)~=-1,
            A=reshape(Factors(FIdx0(c):FIdx1(c)),DimX(c),Fac(c));
            AA=A'*A;
            ssAA=sum(sum(AA.^2));
            ssdiagAA=sum(sum(diag(AA).^2));
            if abs(ssAA-ssdiagAA) > 100*eps;
                Orthomode=0;
            end;
        end;
    end;
end;

if Orthomode==0,
    Zmi=prod(abs(Fac_orig(2:N)));
    Zmj=prod(DimX(2:N));
    Zm=zeros(Zmi,Zmj);
    DimXprodc0 = 1;
    Facprodc0 = 1;
    Zm(1:Facprodc0,1:DimXprodc0)=ones(Facprodc0,DimXprodc0);
    for c=2:N,
        if Fac_orig(c)~=-1,
            A=reshape(Factors(FIdx0(c):FIdx1(c)),DimX(c),Fac(c));
            DimXprodc1 = DimXprodc0*DimX(c);
            Facprodc1 = Facprodc0*Fac(c);
            Zm(1:Facprodc1,1:DimXprodc1)=ckron(A',Zm(1:Facprodc0,1:DimXprodc0));
            DimXprodc0 = DimXprodc1;
            Facprodc0 = Facprodc1;
        end;
    end;
    if Fac_orig(1)~=-1,
        A=reshape(Factors(FIdx0(1):FIdx1(1)),DimX(1),Fac(1));
        Xm=A*G*Zm;
    else
        Xm=G*Zm;
    end;
elseif Orthomode==1,
    CurDimX=DimG;
    Xm=G;
    newi=CurDimX(2);
    newj=prod(CurDimX)/CurDimX(2);
    Xm=reshape(Xm',newi,newj);
    for c=2:N,
        if Fac_orig(c)~=-1,
            A=reshape(Factors(FIdx0(c):FIdx1(c)),DimX(c),Fac(c));
            Xm=A*Xm;
            CurDimX(c)=DimX(c);
        else
            CurDimX(c)=DimX(c);
        end;
        if c~=N,
            newi=CurDimX(c+1);
            newj=prod(CurDimX)/CurDimX(c+1);
        else,
            newi=CurDimX(1);
            newj=prod(CurDimX)/CurDimX(1);
        end;
        Xm=reshape(Xm',newi,newj);
    end;
    if Fac_orig(1)~=-1,
        A=reshape(Factors(FIdx0(1):FIdx1(1)),DimX(1),Fac(1));
        Xm=A*Xm;
    end;
end;

Xm = reshape(Xm,DimX);


function load=pfls(ZtZ,ZtX,dimX,cons,OldLoad,DoWeight,W);

%PFLS
%
% See also:
% 'unimodal' 'monreg' 'fastnnls'
%
%
% Calculate the least squares estimate of
% load in the model X=load*Z' => X' = Z*load'
% given ZtZ and ZtX
% cons defines if an unconstrained solution is estimated (0)
% or an orthogonal (1), a nonnegativity (2), or a unimodality (3)
%
%
% Used by PARAFAC.M

% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $
% Apr 2002 - Fixed error in weighted ls $ rb

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.


if ~DoWeight
    
    if cons==0 % No constr
        %load=((Z'*Z)\Z'*Xinuse)';
        load=(pinv(ZtZ)*ZtX)';
        
    elseif cons==1 % Orthogonal loadings acc. to Harshman & Lundy 94
        load=ZtX'*(ZtX*ZtX')^(-.5);
        
    elseif cons==2 % Nonnegativity constraint
        load=zeros(size(OldLoad));
        for i=1:dimX
            load(i,:)=fastnnls(ZtZ,ZtX(:,i))';
            %       if min(load(i,:))<-eps*1000
            %          load(i,:)=OldLoad(i,:);
            %       end
        end
        
    elseif cons==3 % Unimodality & NNLS
        load=OldLoad;
        F=size(OldLoad,2);
        if F>1
            for i=1:F
                ztz=ZtZ(i,i);
                ztX=ZtX(i,:)-ZtZ(i,[1:i-1 i+1:F])*load(:,[1:i-1 i+1:F])';
                beta=(pinv(ztz)*ztX)';
                load(:,i)=ulsr(beta,1);
            end
        else
            beta=(pinv(ZtZ)*ZtX)';
            load=ulsr(beta,1);
        end
    end
    
elseif DoWeight
    Z=ZtZ;
    X=ZtX;
    if cons==0 % No constr
        load=OldLoad;
        one=ones(1,size(Z,2));
        for i=1:dimX
            ZW=Z.*(W(i,:).^2'*one);
            %load(i,:)=(pinv(Z'*diag(W(i,:))*Z)*(Z'*diag(W(i,:))*X(i,:)'))';
            load(i,:)=(pinv(ZW'*Z)*(ZW'*X(i,:)'))';
        end
        
    elseif cons==2 % Nonnegativity constraint
        load=OldLoad;
        one=ones(1,size(Z,2));
        for i=1:dimX
            ZW=Z.*(W(i,:).^2'*one);
            load(i,:)=fastnnls(ZW'*Z,ZW'*X(i,:)')';
        end
        
    elseif cons==1
        disp(' Weighted orthogonality not implemented yet')
        disp(' Please contact the authors for further information')
        error
        
    elseif cons==3
        disp(' Weighted unimodality not implemented yet')
        disp(' Please contact the authors for further information')
        error
        
    end
    
end


% Check that NNLS and alike do not intermediately produce columns of only zeros
if cons==2|cons==3
    if any(sum(load)==0)  % If a column becomes only zeros the algorithm gets instable, hence the estimate is weighted with the prior estimate. This should circumvent numerical problems during the iterations
        load = .9*load+.1*OldLoad;
    end
end

function pfplot(X,Factors,Weights,Option);
%PFPLOT plot parafac model
%
% See also:
% 'parafac'
%
%
% pfplot(X,Factors,Weights,Option);
% Different aspects for evaluation of the solution.
%
% Option # = 1
% 1	NOT ACCESIBLE
% 2	NOT ACCESIBLE
% 3	DIAGONALITY PLOT
% 4	PLOTS OF RESIDUAL VARIANCE
% 5	PLOTS OF LEVERAGE
% 6	RESIDUALS (STANDARD DEVIATION) VERSUS LEVERAGE
% 7	NORMAL PROBABILITY PLOT
% 8	LOADING PLOT
%
% You HAVE to input all four inputs. If you have no weights, just input [].
% The last input must be an 8-vector with ones if you want the plot and
% zeros else. E.g.
%
% pfplot(X,factors,[],[0 0 1 0 0 0 0 1]);
%
% to have the diagonality and the loading plot
%

% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $
% $ Version 1.03 $ Date 6. October 1999 $ Changed to handle missing values correctly$
% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.



warning off

DimX = size(X);
X = reshape(X,DimX(1),prod(DimX(2:end)));

% Convert to old format
NewLoad = Factors;
ff = [];
for f=1:length(Factors)
    ff=[ff;Factors{f}(:)];
end
Factors = ff;


factors = Factors;
ord=length(DimX);
Fac=length(factors)/sum(DimX);
lidx(1,:)=[1 DimX(1)*Fac];
for i=2:ord
    lidx=[lidx;[lidx(i-1,2)+1 sum(DimX(1:i))*Fac]];
end
if Option(3)==1
    % ESTIMATE DIAGONALITY OF T3-CORE
    diagonality=corcond(reshape(X,DimX),NewLoad,Weights,1);
end
model=nmodel(NewLoad);
model = reshape(model,DimX(1),prod(DimX(2:end)));
if Option(4)==1
    % PLOTS OF RESIDUAL VARIANCE
    figure,eval(['set(gcf,''Name'',''Residual variance'');']);
    aa=ceil(sqrt(ord));bb=ceil(ord/aa);
    for i=1:ord
        r=nshape(reshape(X-model,DimX),i)';
        varian=stdnan(r).^2;
        subplot(aa,bb,i)
        plot(varian)
        if DimX(i)<30
            hold on
            plot(varian,'r+')
        end
        eval(['xlabel(''Mode ', num2str(i),''');']);
        ylabel('Residual variance');
    end
end
if Option(5)==1
    % PLOTS OF LEVERAGE
    figure
    eval(['set(gcf,''Name'',''Leverage'');']);
    aa=ceil(sqrt(ord));
    bb=ceil(ord/aa);
    for i=1:ord
        A=reshape(factors(lidx(i,1):lidx(i,2)),DimX(i),Fac);
        lev=diag(A*pinv(A'*A)*A');
        subplot(aa,bb,i)
        if std(lev)>eps
            plot(lev+100*eps,'+')
            for j=1:DimX(i)
                text(j,lev(j),num2str(j))
            end
        else
            warning('Leverage is constant')
        end
        eval(['xlabel(''Mode ', num2str(i),''');']);
        ylabel('Leverage');
    end
end
if Option(6)==1
    % RESIDUALS (STANDARD DEVIATION) VERSUS LEVERAGE
    figure
    eval(['set(gcf,''Name'',''Residuals vs. Leverages'');']);
    aa=ceil(sqrt(ord));
    bb=ceil(ord/aa);
    for i=1:ord
        subplot(aa,bb,i)
        A=reshape(factors(lidx(i,1):lidx(i,2)),DimX(i),Fac);
        lev=diag(A*pinv(A'*A)*A')'+100*eps;
        r=nshape(reshape(X-model,DimX),i)';
        stand=stdnan(r);
        if std(lev)>eps
            plot(lev,stand,'+')
            for j=1:DimX(i)
                text(lev(j),stand(j),num2str(j))
            end
            eval(['xlabel(''Leverage in mode ', num2str(i),''');']);
            ylabel('Standard deviation');
        else
            warning('Leverage is constant')
        end
    end
end
if Option(7)==1
    % NORMAL PROBABILITY PLOT
    if exist('normplot')
        disp(' ')
        disp(' Normal probability plots are time-consuming')
        disp(' They are made in the statistics toolbox though, so we can''t change that!')
        figure,
        eval(['set(gcf,''Name'',''Normal probability of residuals'');']);
        aa=ceil(sqrt(ord));
        bb=ceil(ord/aa);
        r=nshape(reshape(X-model,DimX),i)';
        r=r(:);
        normplot(r(find(~isnan(r))))
    end
end
if Option(8)==1
    % LOADING PLOT
    if sum(Option)>1
        figure
    end
    eval(['set(gcf,''Name'',''Loadings'');']);
    aa=ceil(sqrt(ord));
    bb=ceil(ord/aa);
    for i=1:ord
        subplot(aa,bb,i)
        A=reshape(factors(lidx(i,1):lidx(i,2)),DimX(i),Fac);
        plot(A)
        eval(['xlabel(''Mode ', num2str(i),''');']);
        ylabel('Loading');
    end
end
drawnow



function id = nident(J,order,mag);

% NIDENT make 'identity' multi-way array
%
% Make 'identity' array of order given by "order" and dimension J.
% id = nident(J,order);
% if extra input vector, mag, is given the j'th superdiagonal will be
% equal to mag(j)

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.


% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $
% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $

if nargin<3
    mag=ones(J,1);
end

id=zeros(J^order,1);
for f=1:J
    idd=f;
    for i=2:order
        idd=idd+J^(i-1)*(f-1);
    end
    id(idd)=mag(f);
end
id  = reshape(id,ones(1,order)*J);

function varargout=nshape(X,f);

%NSHAPE rearrange a multi-way array
%
% [Xf,DimXf] = nshape(X,f);
%
% Refolds an N-way array so that Xf is X with index
% f as row-index, and the remaining in succesive order. For an
% I x J x K x L four-way array this means X1 is I x JKL, X2 is
% J x IKL, X3 is K x IJL, and X4 is L x IJK
%
%
%    K  _______
%      /      /|           1      J     2·J    J·K
%     /______/ |         1  _____________________
%    |      |  |           |      |      |      |
%    |      | /    -->     |      |      |      |        f = (Mode) 1 (same as original array)
% I  |______|/          I  |______|______|______|
%           J
%
%                          1      I     2·I    K·I
%                        1  _____________________
%                          |      |      |      |
%                  -->     |      |      |      |        f = (Mode) 2
%                        J |______|______|______|
%
%
%                          1      I     2·I    I·J
%                        1  _____________________
%                          |      |      |      |
%                  -->     |      |      |      |        f = (Mode) 3
%                        K |______|______|______|
%
%
% f can also indicate the order (meaning the sequence) of the modes
% [Xf,DimXf] = nshape(X,[3 2 1 4]);
% will return Xf as K x JIL
%
% If the last input is not given all rearrangements are given.
% For a fourway array this would read
% [X1,X2,X3,X4]=nshape(X);
%

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.


% $ Version 1.03 $ Date 18. July 1999 $ Not compiled $
% $ Version 1.031 $ Date 18. July 1999 $ Error in help figure and now outputs new DimX $ Not compiled $
% $ Version 2.0 $ Jan 2002 $ Not compiled $ Improved speed and added permute functionality Giorgio Tomasi


ord       = ndims(X);
DimX      = size(X);
varargout = [];
if nargin < 2
    f     = 0;
    do_it = ones(1,nargout);
else
    if length(f) == 1
        do_it = [1:ord] == f;
    end
end
if length(f) == 1
    for i = 1:ord
        if do_it(i)
            varargout{end+1} = reshape(permute(X,[i 1:i-1 i+1:ord]),DimX(i),prod(DimX([1:i-1 i+1:ord])));
        end
    end
    if nargin == 2
        varargout{2} = [DimX(f) DimX([1:f-1 f+1:ord])];
    end
else
    if length(f)==ord
        DimX         = DimX(f);
        varargout{1} = reshape(permute(X,f),DimX(1),prod(DimX(2:end)));
        if nargin == 2
            varargout{2} = DimX;
        end
    else
        error(['f can either be the dimension to be put first or',char(10),...
            'a vector containing the new order of the dimensions']);
    end
end


function [Factors,G,ExplX,Xm]=tucker(X,Fac,Options,ConstrF,ConstrG,Factors,G);
%TUCKER multi-way tucker model
%
% function [Factors,G,ExplX,Xm]=tucker(X,Fac[,Options[,ConstrF,[ConstrG[,Factors[,G]]]]]);
%
% Change: True LS unimodality now supported.
%
% This algorithm requires access to:
% 'fnipals' 'gsm' 'inituck' 'calcore' 'nmodel' 'nonneg' 'setopts' 'misssum'
% 'missmean' 't3core'
%
% See also:
% 'parafac' 'maxvar3' 'maxdia3' 'maxswd3'
%
% ---------------------------------------------------------
%             The general N-way Tucker model
% ---------------------------------------------------------
%
% [Factors,G,ExplX,Xm]=tucker(X,Fac,Options,ConstrF,ConstrG,Factors,G);
% [Factors,G,ExplX,Xm]=tucker(X,Fac);
%
% INPUT
% X        : The multi-way data array.
% Fac      : Row-vector describing the number of factors
%            in each of the N modes. A '-1' (minus one)
%            will tell the algorithm not to estimate factors
%            for this mode, yielding a Tucker2 model.
%            Ex. [3 2 4]
%
% OPTIONAL INPUT
% Options  : See parafac.
% ConstrF  : Constraints that must apply to 'Factors'.
%            Define a row-vector of size N that describes how
%            each mode should be treated.
%            '0' orthogonality (default)
%            '1' non-negativity
%            '2' unconstrained
%            '4' unimodality and non-negativity.
%            E.g.: [0 2 1] yields ortho in first mode, uncon in the second
%            and non-neg in the third mode.
%            Note: The algorithm uses random values if there are no
%            non-negative components in the iteration intermediates. Thus,
%            if non-negativity is applied, the iterations may be
%            non-monotone in minor sequences.
% ConstrG  : Constraints that must apply to 'G'.
%            '[]' or '0' will not constrain the elements of 'G'.
%            To define what core elements should be allowed, give a core that
%            is 1 (one) on all active positions and zero elsewhere - this boolean
%            core array must have the same dimensions as defined by 'Fac'.
%
% OUTPUT
% Factors  : A row-vector containing the solutions.
% G        : Core array that matches the dimensions defined by 'Fac'.
% ExplX    : Fraction of variation (sums of squares explained)
% Xm       : Xhat (the model of X)
%
% This algorithm applies to the general N-way case, so
% the array X can have any number of dimensions. The
% principles of 'projections' and 'systematic unfolding
% methodology (SUM)' are used in this algorithm to provide
% a fast approach - also for larger data arrays. This
% algorithm can handle missing values if denoted
% by NaN's. It can also be used to make TUCKER2/1 models by
% properly setting the elements of 'Fac' to -1.
%
% Note: When estimating a Tucker model on data using non-orthogonal factors,
%       the sum of square of the core may differ between models of the
%       same dataset. This is in order since the factors may
%       thus be correlated. However, the expl. var. should always be the same.
%

% $ Version 2.01  $ May 2011 $ Fixed problem with Tucker2 calculations & missing 'Option'$ CA $ Not compiled $
% $ Version 2.003 $ Jan 2002 $ Fixed problem with length of factors under special conditions $ CA $ Not compiled $
% $ Version 2.002 $ Jan 2002 $ Fixed reshaping of old input G $ RB $ Not compiled $
% $ Version 2.001 $ July 2001 $ Changed problem with checking if Factors exist (should check if it exists in workspace specifically)$ RB $ Not compiled $
% $ Version 2.00  $ May 2001 $ Changed to array notation $ RB $ Not compiled $
% $ Version 1.12  $ Date 14. Nov. 1999 $ Not compiled $
%
%
% Copyright, 1998-, 2011-
% This M-file and the code in it belongs to the holder of the
% copyrights and is made public under the following constraints:
% It must not be changed or modified and code cannot be added.
% The file must be regarded as read-only. Furthermore, the
% code can not be made part of anything but the 'N-way Toolbox'.
% In case of doubt, contact the holder of the copyrights.
%
% Claus A. Andersson
% Chemometrics Group - Life Sciences
% University of Copenhagen
% Rolighedsvej 30, DK-1958 Frederiksberg, Denmark
% E-mail: claus@andersson.dk // ca@life.ku.dk
%

DimX = size(X);
X = reshape(X,DimX(1),prod(DimX(2:end)));
FacNew = Fac;
FacNew(find(FacNew==-1)) = DimX(find(FacNew==-1));

format long
format compact
dbg=0;

if nargin==0,
    help('tucker.m');
    error(['Error calling ''tucker.m''. Since no input arguments were given, the ''help'' command was initiated.'])
    return;
end;
if nargin<2,
    help('tucker.m');
    error(['Error calling ''tucker.m''. At least two (2) input arguments must be given. Read the text above.'])
    return;
end;
if size(Fac,2)==1,
    help('tucker.m');
    error(['Error calling ''tucker.m''. ''Fac'' must be a row-vector.'])
end;

% Initialize system variables
N=size(Fac,2);
Fac_orig=Fac;
finda=find(Fac==-1);
if ~isempty(finda),
    Fac(finda)=zeros(size(finda));
end;
FIdx0=cumsum([1 DimX(1:N-1).*Fac(1:N-1)]);
FIdx1=cumsum([DimX.*Fac]);
pmore=30;
pout=0;
Xm=[];
MissingExist=any(isnan(X(:)));
if MissingExist,
    IdxIsNans=find(isnan(X));
end;
SSX=misssum(misssum(X.^2));

if exist('Options')~=1,
    Options_=[0];
else
    Options_=Options;
end;
Options=[1e-7 0 1 1 1 1000 1 0 1 160];
i=find(Options_);
Options(i)=Options_(i);
if isnan(Options(5)),
    prlvl = 0;
else
    prlvl = 1;
end;
Options12=Options(1);
Options11=Options12*10;
Options21=Options(2);
Options31=Options(3);
Options41=Options(4);
Options51=Options(5);
Options61=Options(6);
Options71=Options(7);
Options81=Options(8);
Options91=Options(9);
Options101=Options(10);

if ~exist('ConstrF'),
    ConstrF=[];
end;
if isempty(ConstrF),
    ConstrF=zeros(size(DimX));
end;
if ConstrF==0 ,
    ConstrF=zeros(size(DimX));
end;

if ~exist('ConstrG')
    ConstrG=[];
end;
if isempty(ConstrG),
    ConstrG=0;
end;

if exist('Factors')~=1,
    Factors=[];
end;

if ~(exist('G')==1),
    G=[];
else
    G=reshape(G,size(G,1),prod(size(G))/size(G,1));
end;

%Give a status/overview
if prlvl>0,
    fprintf('\n\n');
    fprintf('=================   RESUME  &  PARAMETERS   ===================\n');
    fprintf('Array                 : %i-way array with dimensions (%s)\n',N,int2str(DimX));
    if any(Fac==0),
        fprintf('Model                 : (%s) TUCKER2 model\n',int2str(Fac));
    else
        fprintf('Model                 : (%s) TUCKER3 model\n',int2str(Fac));
    end;
end

%Mth initialization
txt1=str2mat('derived by SVD (orthogonality constrained).');
txt1=str2mat(txt1,'derived by NIPALS (orthogonality constrained).');
txt1=str2mat(txt1,'derived by Gram-Schmidt (orthogonality constrained).');
txt1=str2mat(txt1,'This mode is not compressed/calculated, i.e., TUCKER2 model.');
txt1=str2mat(txt1,'derived by non-negativity least squares.');
txt1=str2mat(txt1,'derived by unconstrained simple least squares.');
txt1=str2mat(txt1,'unchanged, left as defined in input ''Factors''.');
txt1=str2mat(txt1,'derived by unimodality constrained regression.');
MethodO=1;
for k=1:N,
    UpdateCore(k)=1;
    if ConstrF(k)==0,
        if Fac(k)>0,
            if 0<DimX(k) & DimX(k)<=180,
                Mth(k)=1;
            end;
            if 180<DimX(k) & DimX(k)<=Inf,
                Mth(k)=3;
            end;
            if Fac(k)<=6 & 180<DimX(k),
                Mth(k)=2;
            end;
        end;
        UpdateWithPinv(k)=1; %Update with the L-LS-P-w/Kron approach
        CalcOrdinar(k)=1;
    end;
    if ConstrF(k)==1,
        Mth(k)=5; %nonneg
        MethodO=2; %use the flexible scheme
        UpdateCore(k)=1; %Update the core in this mode
        CalcOrdinar(k)=1;
    end;
    if ConstrF(k)==2,
        Mth(k)=6; %uncon
        MethodO=2; %use the flexible scheme
        UpdateCore(k)=1; %Update the core in this mode
        CalcOrdinar(k)=1;
    end;
    if ConstrF(k)==3,
        Mth(k)=7; %unchanged
        MethodO=2;
        UpdateCore(k)=1; %Update the core in this mode
        CalcOrdinar(k)=1;
    end;
    if ConstrF(k)==4,
        Mth(k)=8; %unimod
        MethodO=2; %use the flexible scheme
        UpdateCore(k)=1; %Update the core in this mode
        CalcOrdinar(k)=1;
    end;
    if Fac_orig(k)==-1
        Mth(k)=4;
        UpdateCore(k)=0; %Do not update core for this mode
        CalcOrdinar(k)=1;
    end;
    if Options91>=1,
        if prlvl>0,
            if Mth(k)~=4,
                fprintf('Mode %i                : %i factors %s\n',k,Fac(k),txt1(Mth(k),:));
            else
                fprintf('Mode %i                : %s\n',k,txt1(Mth(k),:));
            end;
        end;
    end;
end;

UserFactors=1;
if isempty(Factors),
    UserFactors=0;
else
    ff = [];
    for f=1:length(Factors)
        if ~all(size(Factors{f})==[DimX(f),Fac(f)]), %%
            Factors{f}=rand(DimX(f),Fac(f));%% Added by CA, 27-01-2002
        end;%%
        ff=[ff;Factors{f}(:)];
    end
    OrigFactors=Factors;
    Factors = ff;
end;

usefacinput=0;
if MissingExist,
    if ~UserFactors
        [i j]=find(isnan(X));
        mnx=missmean(X)/3;
        mny=missmean(X')/3;
        n=size(i,1);
        for k=1:n,
            i_=i(k);
            j_=j(k);
            X(i_,j_) = mny(i_) + mnx(j_);
        end;
        mnz=(missmean(mnx)+missmean(mny))/2;
        p=find(isnan(X));
        X(p)=mnz;
    else
        usefacinput=1;
        % Convert to new format
        clear ff,id1 = 0;
        for i = 1:length(DimX)
            id2 = sum(DimX(1:i).*Fac(1:i));ff{i} = reshape(Factors(id1+1:id2),DimX(i),Fac(i));id1 = id2;
        end
        Fact = ff;
        Xm=nmodel(Fact,reshape(G,Fac_orig));
        Xm = reshape(Xm,DimX(1),prod(DimX(2:end)));
        X(IdxIsNans)=Xm(IdxIsNans);
    end;
    SSMisOld=sum(sum( X(IdxIsNans).^2 ));
    SSMis=SSMisOld;
end;

% Initialize the Factors by some method
UserFactors=1;
if isempty(Factors),
    Factors=inituck(reshape(X,DimX),Fac_orig,2,[]);
    
    % Convert to old factors
    ff = [];
    for f=1:length(Factors)
        ff=[ff;Factors{f}(:)];
    end
    Factors = ff;
    UserFactors=0;
end;

% Initialize the core
Core_uncon=0;
Core_nonneg=0;
Core_cmplex=0;
Core_const=0;
G_cons=ConstrG;
if all(ConstrG(:)==1),
    ConstrG=0;
end;
if ConstrG==0,
    % Convert to new format
    clear ff,id1 = 0;
    for i = 1:length(DimX)
        id2 = sum(DimX(1:i).*Fac(1:i));ff{i} = reshape(Factors(id1+1:id2),DimX(i),Fac(i));id1 = id2;
    end
    Fact = ff;
    G=calcore(reshape(X,DimX),Fact,[],1,MissingExist);
    G = reshape(G,size(G,1),prod(size(G))/size(G,1));
    Core_uncon=1;
elseif prod(size(ConstrG)==[Fac(1) prod(Fac(2:N))]),
    tmpM2=1;
    for k=1:N;
        if Mth(k)==4,
            tmpM1=eye(DimX(k));
        else
            tmpM1=reshape(Factors(FIdx0(k):FIdx1(k)),DimX(k),Fac(k));
        end;
        tmpM2=ckron(tmpM2,tmpM1);
    end
    w=ConstrG(:);
    fwz=find(w==1);
    fwnn=find(w==2);
    fwfix=find(w==3);
    G=zeros(prod(Fac),1);
    G(fwz)=tmpM2(:,fwz)\X(:);
    enda=size(Fac,2); %!
    G=reshape(G,Fac(1),prod(Fac(2:enda)));
    Core_cmplex=1;
    MethodO=2;
    UpdateCore=2*ones(1,N);
elseif ConstrG==2,
    % Convert to new format
    clear ff,id1 = 0;
    for i = 1:length(DimX)
        id2 = sum(DimX(1:i).*Fac(1:i));ff{i} = reshape(Factors(id1+1:id2),DimX(i),Fac(i));id1 = id2;
    end
    Fact = ff;
    G=calcore(reshape(X,DimX),Fact,[],1,MissingExist);
    G = reshape(G,size(G,1),prod(size(G))/size(G,1));
    G(find(G<0))=0;
    Core_nonneg=1;
elseif ConstrG==3,
    Core_const=1;
    UpdateCore=0*UpdateCore; %%%Added by CA, 27-01-2002
end;

if prlvl>0
    fprintf('Type of algorithm     : ');
    if MethodO==1,
        fprintf('Orthogonal projections.\n');
    elseif MethodO==2,
        fprintf('Flexible scheme.\n');
    end;
    fprintf('Core                  : ');
    if Core_cmplex==1 & Core_nonneg==0,
        fprintf('Unconstrained composite/sparse core.\n');
    elseif Core_cmplex==1 & Core_nonneg==1,
        fprintf('Non-negativity constrained and composite/sparse core.\n');
    elseif Core_const==1,
        fprintf('Fixed core.\n');
    else
        fprintf('Full unconstrained core.\n');
    end;
end;

if prlvl>0,
    if MissingExist,
        if Options91>=1,
            fprintf('Missing data          : Yes, 2 active loops (expectation maximization).\n');
            fprintf('                        %i values (%.2f%%) out of %i are NaNs/missing.\n',prod(size(IdxIsNans)),100*prod(size(IdxIsNans))/prod(size(X)),prod(size(X)));
            if usefacinput==0,
                fprintf('                        Missing values initialized from column and row means.\n');
            else
                fprintf('                        Missing values initialized from model based on the given input.\n');
            end;
            fprintf('Convergence crit. 1   : %.5g (relative) sum of sq. core elements (corrected for missing values).\n',Options12);
            fprintf('Convergence crit. 2   : %.5g (relative) sum of sq. of pred. missing values.\n',Options11);
            fprintf('Iteration limit       : %i is the maximum number of overall iterations.\n',Options61);
        end;
    else
        if Options91>=1,
            fprintf('Missing data          : No, 1 active loop.\n');
            fprintf('Convergence crit. 1   : %.5g (relative) sum of sq. core elements.\n',Options12);
            fprintf('Iteration limit       : %i is the maximum number of overall iterations.\n',Options61);
        end;
    end;
    fprintf('\n');
    if MissingExist,
        str1=' Iter. 1  |  Corrected sum of   |  Sum of sq. miss.  |   Expl.  ';
        str2='     #    |  sq. core elements  |       values       |  var. [%]';
        fprintf('%s\n',str1);
        fprintf('%s\n\n',str2);
    else
        str1=' Iter. 1  |       Sum of        |   Expl.  ';
        str2='     #    |  sq. core elements  |  var. [%]';
        fprintf('%s\n',str1);
        fprintf('%s\n\n',str2);
    end;
end;
Conv_true=0;
if MethodO==1, %Can use the faster projection technique
    SSGOld=0;
    Converged2=0;
    it1=0;
    itlim1=0;
    t0=clock;
    while ~Converged2,
        Converged1=0;
        while ~Converged1,
            it1=it1+1;
            % Iterate over the modes
            for c=1:N,
                %Compress the data by projections
                if Mth(c)~=4,
                    CurDimX=DimX;
                    RedData=X;
                    for k=1:N;
                        if k~=c,
                            if Mth(k)~=4,
                                kthFactor=reshape(Factors(FIdx0(k):FIdx1(k)),DimX(k),Fac(k));
                                RedData=kthFactor'*RedData;
                                CurDimX(k)=Fac(k);
                            else
                                RedData=RedData;
                            end,
                        end,
                        if k~=N,
                            newi=CurDimX(k+1);
                            newj=prod(CurDimX)/newi;
                        else
                            newi=CurDimX(1);
                            newj=prod(CurDimX)/newi;
                        end;
                        RedData=reshape(RedData',newi,newj);
                    end;
                    %Reshape to the proper unfolding
                    for k=1:(c-1);
                        if k~=c,
                            newi=CurDimX(k+1);
                            newj=prod(CurDimX)/CurDimX(k+1);
                        else,
                            newi=CurDimX(1);
                            newj=prod(CurDimX)/CurDimX(1);
                        end;
                        RedData=reshape(RedData',newi,newj);
                    end;
                    %Find a basis in the projected space
                    %...using the robust SVD
                    if Mth(c)==1,
                        if MissingExist,
                            [U S V]=svd(RedData',0);
                            cthFactor=V(:,1:Fac(c));
                        else
                            [U S V]=svd(RedData',0);
                            cthFactor=V(:,1:min(Fac(c),size(V,2)));
                            if size(cthFactor,2)<Fac(c)
                                cthFactor = [cthFactor rand(size(cthFactor,1),Fac(c)-size(cthFactor,2))];
                            end
                        end;
                    end;
                    %...using the fast NIPALS
                    if Mth(c)==2,
                        if MissingExist,
                            [cthFactor]=fnipals(RedData,Fac(c),reshape(Factors(FIdx0(c):FIdx1(c)),DimX(c),Fac(c)));
                        else
                            [cthFactor]=fnipals(RedData,Fac(c),reshape(Factors(FIdx0(c):FIdx1(c)),DimX(c),Fac(c)));
                        end;
                    end;
                    %...using simplified continuous Gram-Schmidt orthogonalization
                    if Mth(c)==3,
                        if MissingExist,
                            TempMat=RedData*RedData';
                            cthFactor=reshape(Factors(FIdx0(c):FIdx1(c)),DimX(c),Fac(c));
                            for i=1:2,
                                [cthFactor]=gsm(TempMat*cthFactor);
                            end;
                        else
                            TempMat=RedData*RedData';
                            cthFactor=reshape(Factors(FIdx0(c):FIdx1(c)),DimX(c),Fac(c));
                            for i=1:2,
                                [cthFactor]=gsm(TempMat*cthFactor);
                            end;
                        end;
                    end;
                    %...this is void (no compression for this mode)
                    if Mth(c)==4,
                    end;
                    %...this is void (Keep factors unchanged)
                    if Mth(c)==7,
                    end;
                    %Update the 'Factors' with the current estimates
                    if Mth(c)~=4 & Mth(c)~=7
                        Factors(FIdx0(c):FIdx1(c))=cthFactor(:)';
                    end;
                end;
            end;
            
            if ~Core_const & Core_uncon==1 & Core_nonneg==0,
                % Convert to new format
                clear ff,id1 = 0;
                for i = 1:length(DimX)
                    id2 = sum(DimX(1:i).*Fac(1:i));ff{i} = reshape(Factors(id1+1:id2),DimX(i),Fac(i));id1 = id2;
                end
                Fact = ff;
                G=calcore(reshape(X,DimX),Fact,[],1,MissingExist);
                G = reshape(G,size(G,1),prod(size(G))/size(G,1));
            elseif Core_nonneg==1,
                g=t3core(reshape(X,DimX),Fac,Factors(:),0,1);
                G=reshape(g,Fac(1),prod(Fac(2:N)));
            else
                tmpM2=1;
                for k=1:N;
                    if Mth(k)==4,
                        tmpM1=eye(DimX(k));
                    else
                        tmpM1=reshape(Factors(FIdx0(k):FIdx1(k)),DimX(k),Fac(k));
                    end;
                    tmpM2=ckron(tmpM2,tmpM1);
                end
                G=G(:);
                G(fwz)=tmpM2(:,fwz)\X(:);
                enda=size(Fac,2);
                G=reshape(G,Fac(1),prod(Fac(2:enda)));
            end;
            
            SSG=sum(sum(G.^2));
            if MissingExist,
                SSG=SSG-SSMis;
            end;
            if abs(SSG-SSGOld)<Options12*SSGOld,
                Converged1=1;
            end;
            if it1>=Options61,
                itlim1=1;
                Converged1=1;
                Converged2=1;
            end;
            SSGOld=SSG;
            js=0;
            %Save on time count
            if Options101>0 & (etime(clock,t0)>Options101),
                save('temp.mat','Factors','G','DimX','Fac');
                t0=clock;
                js=1;
            end;
            %Save on iteration count
            %if (Options101<0) & (mod(it1,abs(Options101))==0),
            keval=it1/Options51;
            if (Options101<0) & ( abs( keval - floor(keval) ) <=eps),
                save('temp.mat','Factors','G','DimX','Fac');
                js=1;
            end;
            %if mod(it1,Options51)==0 | it1==1  | js==1,
            keval=it1/Options51;
            if (abs( keval - floor(keval) ) <=eps) | it1==1  | js==1, %Matlab 4.2 comp.
                % Convert to new format
                clear ff,id1 = 0;
                for i = 1:length(DimX)
                    if Fac(i)
                        id2 = sum(DimX(1:i).*Fac(1:i).*(Fac(1:i)~=0));
                        ff{i} = reshape(Factors(id1+1:id2),DimX(i),Fac(i));id1 = id2;
                    else
                        ff{i}=[];
                    end
                end
                Fact = ff;
                Xm=nmodel(Fact,reshape(G,FacNew));
                Xm = reshape(Xm,DimX(1),prod(DimX(2:end)));
                
                if MissingExist,
                    X(IdxIsNans)=Xm(IdxIsNans);
                    SSMis=sum(sum( Xm(IdxIsNans).^2 ));
                    if abs(SSMis-SSMisOld)<Options11*SSMisOld,
                        Converged2=1;
                    end;
                    SSMisOld=SSMis;
                else
                    Converged2=1;
                end;
                ExplX=100*(1-sum(sum((X-Xm).^2))/SSX);
                pout=pout+1;
                if pout>pmore,
                    if prlvl > 0,
                        fprintf('%s\n',str1);
                        fprintf('%s\n',str2);
                    end;
                    pout=0;
                end;
                
                if prlvl>0,
                    if MissingExist,
                        fprintf(' %6i       %14.3f     %14.3f         %8.4f',it1,SSG,SSMis,ExplX);
                    else
                        fprintf(' %6i        %14.3f      %8.4f',it1,SSG,ExplX);
                    end;
                    if js,
                        fprintf(' - saved to ''temp.mat'' \n')
                    else
                        fprintf('\n')
                    end;
                end;
            end;
        end; %Inner loop
    end; %Outer loop
    if prlvl>0,
        if itlim1==0,
            fprintf('   Stopped. Convergence criteria reached.\n');
        else
            fprintf('   Stopped. Iteration limits reached in model and expectation loops.\n');
        end;
        if MissingExist,
            fprintf(' %6i       %14.3f     %14.3f         %8.4f',it1,SSG,SSMis,ExplX);
        else
            fprintf(' %6i        %14.3f      %8.4f',it1,SSG,ExplX);
        end;
    end;
    if Options101~=0,
        save('temp.mat','Factors','G','DimX','Fac');
        if prlvl>0,
            fprintf(' - saved to ''temp.mat'' \n')
        end;
    else
        if prlvl>0,
            fprintf('\n')
        end;
    end;
    
elseif MethodO==2, %Must use slower but more general schemes
    
    SSGOld=0;
    OldFactors=0*Factors;
    Converged2=0;
    it1=0;
    t0=clock;
    itlim1=0;
    while ~Converged2,
        Converged1=0;
        while ~Converged1,
            it1=it1+1;
            %Iterate over the modes
            for c=1:N,
                faclist1=[N:-1:1 N:-1:1];
                faclist=faclist1(N-c+2:N-c+2+(N-2));
                tmpM2=1;
                tmpM2Pinv=1;
                for k=1:N-1;
                    if Mth(faclist(k))==4,
                        tmpM1=eye(DimX(faclist(k)));
                    else
                        tmpM1=reshape(Factors(FIdx0(faclist(k)):FIdx1(faclist(k))),DimX(faclist(k)),Fac(faclist(k)));
                    end;
                    if CalcOrdinar(c)==1,
                        tmpM2=ckron(tmpM2,tmpM1');
                    end
                end
                %Estimate Factors for the cth way
                %...this is void (no compression for this mode)
                if any(Mth(c)==[1 2 3]),
                    tmpM4=G*tmpM2;
                    if MissingExist,
                        %cthFactor=(X*tmpM4')/((tmpM4*X'*X*tmpM4')^(1/2));
                        [U S V]=svd(tmpM4*X'*X*tmpM4',0);
                        mS=min(size(S));
                        Sm=S;
                        Sm(1:mS,1:mS)=diag(1./sqrt(diag(S(1:mS,1:mS))));
                        pinvm=U*Sm*V';
                        cthFactor=X*tmpM4'*pinvm;
                    else
                        %cthFactor=(X*tmpM4')/((tmpM4*X'*X*tmpM4')^(1/2));
                        [U S V]=svd(tmpM4*X'*X*tmpM4',0);
                        mS=min(size(S));
                        Sm=S;
                        Sm(1:mS,1:mS)=diag(1./sqrt(diag(S(1:mS,1:mS))));
                        pinvm=U*Sm*V';
                        cthFactor=X*tmpM4'*pinvm;
                    end;
                end;
                if Mth(c)==4,
                end;
                if Mth(c)==5, %Nonneg
                    tmpM3=reshape(Factors(FIdx0(c):FIdx1(c)),DimX(c),Fac(c));
                    if it1==1,
                        %tmpM3(find(tmpM3<0))=0;
                        i__=find(tmpM3<0);
                        if ~isempty(i__),
                            tmpM3(i__)=zeros(1,length(i__));
                        end;
                    end;
                    if CalcOrdinar(c) == 1,
                        tmpM4=G*tmpM2;
                        if dbg, ss1=sum(sum( (X-tmpM3*tmpM4).^2 ));end;
                        [cthFactor,SS]=nonneg(X,tmpM4,tmpM3);
                        if dbg, ss2=sum(sum( (X-cthFactor*tmpM4).^2 ));
                            fprintf('Nonng report (Ordi) %15.8d  %15.8d\n',ss1,ss2); end;
                    end;
                end;
                if Mth(c)==6, %Uncon
                    cthFactor=reshape(Factors(FIdx0(c):FIdx1(c)),DimX(c),Fac(c));
                    if CalcOrdinar(c) == 1,
                        M_=G*tmpM2;
                        if dbg, ss1=sum(sum( (X-cthFactor*G*tmpM2).^2 ));end;
                        cthFactor=X/M_;
                        if dbg, ss2=sum(sum( (X-cthFactor*G*tmpM2).^2 ));
                            fprintf('Uncon report (Ordi) %15.8d  %15.8d\n',ss1,ss2);end;
                    end;
                end;
                if Mth(c)==7, %Not updated
                end;
                if Mth(c)==8, %Unimodality
                    tmpM3=reshape(Factors(FIdx0(c):FIdx1(c)),DimX(c),Fac(c));
                    if it1==1,
                        %tmpM3(find(tmpM3<0))=0;
                        i__=find(tmpM3<0);
                        if ~isempty(i__),
                            tmpM3(i__)=zeros(1,length(i__));
                        end;
                    end;
                    if CalcOrdinar(c) == 1,
                        tmpM4=G*tmpM2;
                        if dbg, ss1=sum(sum( (X-tmpM3*tmpM4).^2 ));end;
                        cthFactor = unimodal(tmpM4',X');
                        if dbg, ss2=sum(sum( (X-cthFactor*tmpM4).^2 ));fprintf('Unimod report (Ordi) %15.8d  %15.8d\n',ss1,ss2); end;
                    end;
                end;
                
                
                %Update 'Factors' with the current estimates
                if Mth(c)~=4 & Mth(c)~=7
                    cn=sum(cthFactor.^2);
                    for i=1:Fac(c);
                        if it1<=10 & cn(i)<eps,
                            cthFactor(:,i)=rand(size(cthFactor,1),1);
                        end;
                    end;
                    if Core_const,
                        if c>1
                            cthFactor=normit(cthFactor);
                        end;
                    else
                        cthFactor=normit(cthFactor);
                    end;
                    Factors(FIdx0(c):FIdx1(c))=cthFactor(:)';
                end;
                
                %Estimate the new core after each SUBiteration
                CoreupdatedInner=0;
                if ~Core_const
                    if UpdateCore(c)==1,
                        tmpM3=reshape(Factors(FIdx0(c):FIdx1(c)),DimX(c),Fac(c));
                        if CalcOrdinar(c) == 1,
                            if dbg, ss1=sum(sum( (X-tmpM3*G*tmpM2).^2 ));end;
                            G=tmpM3\(X/tmpM2);
                            if dbg, ss2=sum(sum( (X-tmpM3*G*tmpM2).^2 ));
                                fprintf('Core report  (Ordi) %15.8d  %15.8d\n',ss1,ss2);end;
                        end;
                        CoreupdatedInner=1;
                    end;
                    if UpdateCore(c)==2,
                        G=zeros(Fac(c),prod(Fac(nsetdiff(1:N,c))));
                        tmpM2=1;
                        for k=[(c-1):-1:1 N:-1:c];
                            if Mth(k)==4,
                                tmpM1=eye(DimX(k));
                            else
                                tmpM1=reshape(Factors(FIdx0(k):FIdx1(k)),DimX(k),Fac(k));
                            end;
                            tmpM2=ckron(tmpM2,tmpM1);
                        end
                        G=G(:);
                        fwz=find(G_cons==1);
                        G(fwz)=tmpM2(:,fwz)\X(:);
                        G=reshape(G,Fac(c),prod(Fac(nsetdiff(1:N,c))));
                        CoreupdatedInner=1;
                    end;
                end;
                
                %Reshape to the next unfolding
                if c~=N,
                    newix=DimX(c+1);
                    newjx=prod(DimX)/DimX(c+1);
                    newig=FacNew(c+1);
                    newjg=prod(FacNew)/FacNew(c+1);
                else
                    newix=DimX(1);
                    newjx=prod(DimX)/DimX(1);
                    newig=FacNew(1);
                    newjg=prod(FacNew)/FacNew(1);
                end;
                X=reshape(X',newix,newjx);
                G=reshape(G',newig,newjg);
                if length(G_cons(:))>1,
                    G_cons=reshape(G_cons',newig,newjg);
                end;
            end;
            
            %Estimate the new core after each MAIN iteration
            if CoreupdatedInner==0,
                if ~Core_const & Core_cmplex==0,
                    if Core_nonneg==1,
                        g=t3core(reshape(X,DimX),Fac,Factors,0,1);
                        G=reshape(g,Fac(1),prod(Fac(2:N)));
                    elseif Core_nonneg==0,
                        % Convert to new format
                        clear ff,id1 = 0;
                        for i = 1:length(DimX)
                            id2 = sum(DimX(1:i).*Fac(1:i));ff{i} = reshape(Factors(id1+1:id2),DimX(i),Fac(i));id1 = id2;
                        end
                        Fact = ff;
                        G=calcore(reshape(X,DimX),Fact,[],0,MissingExist);
                        G = reshape(G,size(G,1),prod(size(G))/size(G,1));
                    end
                end;
            end;
            CoreupdatedInner=0;
            
            if ~Core_const
                SSG=sum(sum(G.^2));
                if MissingExist,
                    SSG=SSG-SSMis;
                end;
                if abs(SSG-SSGOld)<Options12*SSGOld,
                    Converged1=1;
                end;
                if it1>=Options61,
                    Converged1=1;
                    Converged2=1;
                    itlim1=1;
                end;
                SSGOld=SSG;
            else
                SSF=sum(sum((Factors-OldFactors).^2));
                if SSF<Options12*sum(sum(Factors.^2));
                    Converged1=1;
                end;
                if it1>=Options61,
                    Converged1=1;
                    Converged2=1;
                    itlim1=1;
                end;
                OldFactors=Factors;
                SSG=SSF;
            end;
            
            js=0;
            if Options101>0 & (etime(clock,t0)>Options101),
                save('temp.mat','Factors','G','DimX','Fac');
                t0=clock;
                js=1;
            end;
            %if Options101<0 & mod(it1,abs(Options101))==0,
            keval=it1/Options101;
            if (Options101<0) & (abs( keval - floor(keval) ) <=eps), %Matlab 4.2 comp.
                save('temp.mat','Factors','G','DimX','Fac');
                js=1;
            end;
            %if mod(it1,Options51)==0 | it1==1 | js==1,
            keval=it1/Options51;
            if (abs( keval - floor(keval) ) <=eps) | it1==1  | js==1, %Matlab 4.2 comp.
                % Convert to new format
                clear ff
                id1 = 0;
                for i = 1:length(DimX)
                    if Fac(i)>0
                        id2 = sum(DimX(1:i).*Fac(1:i));
                        ff{i} = reshape(Factors(id1+1:id2),DimX(i),FacNew(i));
                        id1 = id2;
                    end
                end
                Fact = ff;
                Xm = nmodel(Fact,reshape(G,FacNew));
                Xm = reshape(Xm,DimX(1),prod(DimX(2:end)));
                if MissingExist,
                    X(IdxIsNans)=Xm(IdxIsNans);
                    SSMis=sum(sum( Xm(IdxIsNans).^2 ));
                    if abs(SSMis-SSMisOld)<Options11*SSMisOld,
                        Converged2=1;
                    end;
                    SSMisOld=SSMis;
                else
                    Converged2=1;
                end;
                ExplX=100*(1 - sum(sum((X-Xm).^2))/SSX );
                pout=pout+1;
                if pout>pmore,
                    if prlvl>0,
                        fprintf('%s\n',str1);
                        fprintf('%s\n',str2);
                    end;
                    pout=0;
                end;
                if prlvl>0,
                    if MissingExist,
                        fprintf(' %6i       %14.3f     %14.3f         %8.4f',it1,SSG,SSMis,ExplX);
                    else
                        fprintf(' %6i        %14.3f      %8.4f',it1,SSG,ExplX);
                    end;
                    if js,
                        fprintf(' - saved to ''temp.mat'' \n')
                    else
                        fprintf('\n')
                    end;
                end;
            end;
        end; %Inner loop
    end; %Outer loop
    if prlvl>0,
        if itlim1==0,
            fprintf('   Stopped. Convergence criteria reached.\n');
        else
            fprintf('   Stopped. Iteration limits reached in model and expectation loops.\n');
        end;
        if MissingExist,
            fprintf(' %6i       %14.3f     %14.3f         %8.4f',it1,SSG,SSMis,ExplX);
        else
            fprintf(' %6i        %14.3f      %8.4f',it1,SSG,ExplX);
        end;
    end;
    if Options101~=0,
        save('temp.mat','Factors','G','DimX','Fac');
        if prlvl>0,
            fprintf(' - saved to ''temp.mat'' \n')
        end;
    else
        if prlvl>0,
            fprintf('\n')
        end;
    end;
end; %Outer loop

if MissingExist,
    Xf=Xm(IdxIsNans);
end;

Factors=Factors';
format
Xm = reshape(Xm,DimX);
G = reshape(G,FacNew);

% Convert to new format
clear ff
id1 = 0;
for i = 1:length(DimX)
    if Fac(i)
        id2 = sum(DimX(1:i).*Fac(1:i));
        ff{i} = reshape(Factors(id1+1:id2),DimX(i),Fac(i));
        id1 = id2;
    else
        ff{i}=[];
    end
end

Factors = ff;



function [T,P]=fnipals(X,w,T)

%FNIPALS nipals algorithm for PCA
%
% function [T,P]=fnipals(X,w,T)
%
% 'fnipals.m'
%
% This algorithm requires the presence of:
% 'missmean.m'
%
% ----------------------------------------------------
%        Find eigenvectors according to NIPALS
% ----------------------------------------------------
%
% [T,P]=fnipals(X,w,T);
% [T,P]=fnipals(X,w);
%
% T is found so that X = T*P', s.t ||T||=1 and T'T=I
%
% X        : The matrix to be decomposed.
% w        : Number of factors to extract.
%            If w is high (perhaps>20) consider using SVD.
% T        : Initial guess of the solution, optional.
%            If T is not specified, a little time will
%            be used on finding orthogonal random
%            starting values.
%
% You may want to calculate P afterwards by typing 'P=X*T'.
% Note that the T returned is orthonormal.
% Calculation of P is left of this implementation to save FLOP's.
% It handles missing values NaNs (very dispersed, less than 15%)
% If the problem is small enough you would prefer the SVD rather
% than NIPALS for finding T. NIPALS may be inaccurate when
% extracting too many factors, i.e., many more than the rank
% of X.

%scalar ConvLim WarnLim ItMax a b i

% $ Version 1.01 $ Date 18. June 1998 $ Not compiled $

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.


ConvLim=1e-12;
WarnLim=1e-4;
ConvLimMiss=100*ConvLim;
ItMax=100;

filename='fnipals.m';

[a b]=size(X);

if (w>a | w>b) | w<1,
    help(filename);
    error(['Error in ' filename ': Number of factors to extract is invalid!'])
end;

np=isnan(X);
MissingExist=any(np);

if ~exist('T'),
    T=orth(randn(a,w));
end;

if exist('P'),
    P=[];
end;

if ~MissingExist
    if (size(T) == [a w]),
        if a>b,
            P=X'*T;
            l2=Inf;
            Z=X'*X;
            for i=1:w,
                p=P(:,i);
                d=1;
                it=0;
                while (d>ConvLim) & (it<ItMax),
                    it=it+1;
                    p=Z*p;
                    l1=sqrt(p'*p);
                    p=p/l1;
                    d=(l1-l2)^2;
                    l2=l1;
                end;
                P(:,i)=sqrt(l1)*p;
                Z=Z-P(:,i)*P(:,i)';
                WarnLim=sqrt(l1)/1000;
                if it>=ItMax & d>WarnLim,
                    disp('FNIPALS, High-X: Iterated up to the ItMax limit!')
                    disp('FNIPALS, High-X: The solution has not converged!')
                end;
            end;
            T=X*P;
        else
            P=[];
            l2=Inf;
            Z=X*X';
            for i=1:w,
                t=T(:,i);
                d=1;
                it=0;
                while (d>ConvLim) & (it<ItMax),
                    it=it+1;
                    t=Z*t;
                    l1=sqrt(t'*t);
                    t=t/l1;
                    d=(l1-l2).^2;
                    l2=l1;
                end;
                T(:,i)=sqrt(l1)*t;
                Z=Z-T(:,i)*T(:,i)';
                WarnLim=sqrt(l1)/1000;
                if it>=ItMax & d>WarnLim,
                    disp('FNIPALS, Wide-X: Iterated up to the ItMax limit!')
                    disp('FNIPALS, Wide-X: The solution has not converged!')
                end;
            end;
        end;
        T=gsm(T);
    else
        error(['Error in ' filename ': Number of factors to extract is invalid!'])
    end;
else
    MissIdx=find(np);
    [i j]=find(np);
    mnx=missmean(X)/2;
    mny=missmean(X')/2;
    n=size(i,1);
    for k=1:n,
        i_i=i(k);
        j_j=j(k);
        X(i_i,j_j) = mny(i_i) + mnx(j_j);
    end;
    mnz=(missmean(mnx)+missmean(mny))/2;
    
    ssmisold=sum(sum( X(MissIdx).^2 ));
    sstotold=sum(sum( X.^2 ));
    ssrealold=sstotold-ssmisold;
    iterate=1;
    while iterate
        
        if (size(T) == [a w]),
            if a>b,
                P=X'*T;
                l2=Inf;
                Z=X'*X;
                for i=1:w,
                    p=P(:,i);
                    d=1;
                    it=0;
                    while (d>ConvLim) & (it<ItMax),
                        it=it+1;
                        p=Z*p;
                        l1=sqrt(p'*p);
                        p=p/l1;
                        d=(l1-l2)^2;
                        l2=l1;
                    end;
                    P(:,i)=sqrt(l1)*p;
                    Z=Z-P(:,i)*P(:,i)';
                    WarnLim=sqrt(l1)/1000;
                    if it>=ItMax & d>WarnLim,
                        disp('FNIPALS, High-X: Iterated up to the ItMax limit!')
                        disp('FNIPALS, High-X: The solution has not converged!')
                    end;
                end;
                T=X*P;
            else
                P=[];
                l2=Inf;
                Z=X*X';
                for i=1:w,
                    t=T(:,i);
                    d=1;
                    it=0;
                    while (d>ConvLim) & (it<ItMax),
                        it=it+1;
                        t=Z*t;
                        l1=sqrt(t'*t);
                        t=t/l1;
                        d=(l1-l2).^2;
                        l2=l1;
                    end;
                    T(:,i)=sqrt(l1)*t;
                    Z=Z-T(:,i)*T(:,i)';
                    WarnLim=sqrt(l1)/1000;
                    if it>=ItMax & d>WarnLim,
                        disp('FNIPALS, Wide-X: Iterated up to the ItMax limit!')
                        disp('FNIPALS, Wide-X: The solution has not converged!')
                    end;
                end;
            end;
            T=gsm(T);
        else
            error(['Error in ' filename ': Number of factors to extract is invalid!'])
        end;
        
        P=X'*T;
        Xm=T*P';
        X(MissIdx)=Xm(MissIdx);
        ssmis=sum(sum( Xm(MissIdx).^2 ));
        sstot=sum(sum( X.^2 ));
        ssreal=sstot-ssmis;
        if abs(ssreal-ssrealold)<ConvLim*ssrealold & abs(ssmis-ssmisold)<ConvLimMiss*ssmisold,
            iterate=0;
        end;
        ssrealold=ssreal;
        ssmisold=ssmis;
    end;
end;
T=gsm(T);


function [E]=GSM(V);
%GSM orthogonalization
%
% [E]=GSM(V);
% GS   Gram-Schmidt Method for orthogonalisation
%      An orthonormal basis spanning the columns of V is returned in E.
%
%      This algorithm does not use pivoting or any other
%      stabilization scheme. For a completely safe orthogonalization
%      you should use 'ORTH()' though is may take triple the time.
%      'GSM()' is optimized for speed and requies only minimum storage
%      during iterations. No check of rank is performed on V!
%

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.


[m n]=size(V);

%Allocate space for the basis
E=zeros(m,n);

%The first basis vector is taken directly from V
s=sqrt(sum(V(:,1).^2));
E(:,1)=V(:,1)/s;

%Find the other basis vectors as orthogonals to
%the already determined basis by projection
for k=2:n,
    f=V(:,k)-E(:,1:(k-1))*(E(:,1:(k-1))'*V(:,k));
    s=sqrt(sum(f.^2));
    if s<eps,
        E(:,k)=0*f;   %set to zeros
    else
        E(:,k)=f/s;   %normalize
    end;
end;

function [x,w] = fastnnls(XtX,Xty,tol)
%NNLS	Non-negative least-squares.
%	b = fastnnls(XtX,Xty) returns the vector b that solves X*b = y
%	in a least squares sense, subject to b >= 0, given the inputs
%       XtX = X'*X and Xty = X'*y.
%
%	A default tolerance of TOL = MAX(SIZE(X)) * NORM(X,1) * EPS
%	is used for deciding when elements of b are less than zero.
%	This can be overridden with b = fastnnls(X,y,TOL).
%
%	[b,w] = fastnnls(XtX,Xty) also returns dual vector w where
%	w(i) < 0 where b(i) = 0 and w(i) = 0 where b(i) > 0.
%
%
%	L. Shure 5-8-87 Copyright (c) 1984-94 by The MathWorks, Inc.
%
%  Revised by:
%	Copyright
%	Rasmus Bro 1995
%	Denmark
%	E-mail rb@kvl.dk
%  According to Bro & de Jong, J. Chemom, 1997, 11, 393-401

% initialize variables


if nargin < 3
    tol = 10*eps*norm(XtX,1)*max(size(XtX));
end
[m,n] = size(XtX);
P = zeros(1,n);
Z = 1:n;
x = P';
ZZ=Z;
w = Xty-XtX*x;

% set up iteration criterion
iter = 0;
itmax = 30*n;

% outer loop to put variables into set to hold positive coefficients
while any(Z) & any(w(ZZ) > tol)
    [wt,t] = max(w(ZZ));
    t = ZZ(t);
    P(1,t) = t;
    Z(t) = 0;
    PP = find(P);
    ZZ = find(Z);
    nzz = size(ZZ);
    z(PP')=(Xty(PP)'/XtX(PP,PP)');
    z(ZZ) = zeros(nzz(2),nzz(1))';
    z=z(:);
    % inner loop to remove elements from the positive set which no longer belong
    
    while any((z(PP) <= tol)) & iter < itmax
        
        iter = iter + 1;
        QQ = find((z <= tol) & P');
        alpha = min(x(QQ)./(x(QQ) - z(QQ)));
        x = x + alpha*(z - x);
        ij = find(abs(x) < tol & P' ~= 0);
        Z(ij)=ij';
        P(ij)=zeros(1,max(size(ij)));
        PP = find(P);
        ZZ = find(Z);
        nzz = size(ZZ);
        z(PP)=(Xty(PP)'/XtX(PP,PP)');
        z(ZZ) = zeros(nzz(2),nzz(1));
        z=z(:);
    end
    x = z;
    w = Xty-XtX*x;
end

x=x(:);


function G=t3core(X,Load,Weights,NonNeg);
%T3CORE calculate Tucker core
%
% G=T3core(X,Load,Weights,NonNeg);
% Calculate a Tucker3 core given X, the loadings, Load
% in vectorized format and optionally Weights. Missing NaN, NonNeg = 1 => nonnegativity


% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.


% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $




for i = 1:length(Load)
    Fac(i) = size(Load{i},2);
end

DimX = size(X);
X = reshape(X,DimX(1),prod(DimX(2:end)));

% Convert to old format
ff = [];
for f=1:length(Load)
    ff=[ff;Load{f}(:)];
end
Load = ff;

ord=length(DimX);

if length(Fac)==1
    Fac=ones(1,ord)*Fac;
end

l_idx=0;
for i=1:ord
    l_idx=[l_idx sum(DimX(1:i).*Fac(1:i))];
end

if exist('Weights')~=1,
    Weights=0;
elseif isempty(Weights)
    Weights = 0;
end

if exist('NonNeg')~=1,
    NonNeg=0;
end


if any(isnan(X(:)))
    id=find(isnan(X));
    M=ones(size(X));
    M(id)=zeros(size(id));
    X(id)=1000000*ones(size(id));
    if prod(size(Weights))==prod(DimX)
        %Modified by GT
        Weights = reshape(Weights,DimX(1),prod(DimX(2:end)));
        %End GT
        Weights = Weights.*M;
    else
        Weights=M;
    end
end


if prod(size(Weights))==prod(DimX) % Use weighted approach
    LL=reshape(Load(l_idx(1)+1:l_idx(2)),DimX(1),Fac(1));
    xtz=zeros(1,prod(Fac));
    ztz=zeros(prod(Fac),prod(Fac));
    for i=1:DimX(1)
        L1=reshape(Load(l_idx(ord)+1:l_idx(ord+1)),DimX(ord),Fac(ord));
        L2=reshape(Load(l_idx(ord-1)+1:l_idx(ord)),DimX(ord-1),Fac(ord));
        Z=kron(L1,L2);
        for ii=ord-2:-1:2
            L=reshape(Load(l_idx(ii)+1:l_idx(ii+1)),DimX(ii),Fac(ii));
            Z=kron(Z,L);
        end
        Z=kron(Z,LL(i,:));
        ztz=ztz+(Z.*(Weights(i,:)'*ones(1,size(Z,2)) ))'*Z;
        xtz=xtz+(X(i,:).*Weights(i,:))*Z;
    end
    if NonNeg==1;
        G=fastnnls(ztz,xtz');
    else
        G=pinv(ztz)*xtz';
    end
    
else % No weighting
    
    ztz=zeros(prod(Fac),prod(Fac));
    
    L1=reshape(Load(l_idx(ord)+1:l_idx(ord+1)),DimX(ord),Fac(ord));
    L1tL1=L1'*L1;
    L2=reshape(Load(l_idx(ord-1)+1:l_idx(ord)),DimX(ord-1),Fac(ord-1));
    L2tL2=L2'*L2;
    ztz=kron(L1tL1,L2tL2);
    for o=ord-2:-1:1,
        L=reshape(Load(l_idx(o)+1:l_idx(o+1)),DimX(o),Fac(o));
        LtL=L'*L;
        ztz=kron(ztz,LtL);
    end
    
    xtz=zeros(prod(Fac),1);
    F=ones(ord,1);
    F(1)=0;
    for f=1:prod(Fac)
        F(1)=F(1)+1;
        for ff=1:ord-1
            if F(ff)==Fac(ff)+1;
                F(ff+1)=F(ff+1)+1;
                F(ff)=1;
            end
            % F runs through all combinations of factors
        end
        L=reshape(Load(l_idx(1)+1:l_idx(2),:),DimX(1),Fac(1));
        cc=L(:,F(1))'*X;
        for j=ord:-1:3,
            ccc=zeros(prod(DimX(2:j-1)),DimX(j));
            ccc(:)=cc;
            cc=ccc';
            L=reshape(Load(l_idx(j)+1:l_idx(j+1),:),DimX(j),Fac(j));
            cc=L(:,F(j))'*cc;
        end
        L=reshape(Load(l_idx(2)+1:l_idx(2+1),:),DimX(2),Fac(2));
        cc=L(:,F(2))'*cc';
        xtz(f)=cc;
    end
    if NonNeg==1,
        G=fastnnls(ztz,xtz);
    else
        G=pinv(ztz)*xtz;
    end
end

G = reshape(G,Fac);

function [A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,X,Y,Z]=fac2let(Factors,DimX);

%FAC2LET Convert 'Factors' to component matrices
%
%
% [A,B,C]=fac2let(Factors);
% [A,B,C,D]=fac2let(Factors);
% [A,B,C,D,E]=fac2let(Factors);
%             .....
% [A,B,C,...,Z]=fac2let(Factors);
%
% This algorithm applies to the N-way case (2<N<25).
%

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.


% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $
% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $

Txt='ABCDEFGHIJKLMNOPQRSTUVXYZ';

if nargin==2
    order = length(DimX);
    F     = prod(length(Factors))/sum(DimX);
    for i=1:order
        start = sum(DimX(1:i-1))*F+1;
        endd = sum(DimX(1:i))*F;
        eval([Txt(i) ,'= reshape(Factors(',num2str(start),':',num2str(endd),'),',num2str(DimX(i)),',',num2str(F),');']);
    end
else
    for i = 1:length(Factors)
        eval([Txt(i) ,'= Factors{i};']);
    end
end

function [A,B,C]=gram(X1,X2,F);

%GRAM generalized rank annihilation method
%
% [A,B,C]=gram(X1,X2,F);
%
% cGRAM - Complex Generalized Rank Annihilation Method
% Fits the PARAFAC model directly for the case of a
% three-way array with only two frontal slabs.
% For noise-free trilinear data the algorithm is exact.
% If input is not complex, similarity transformations
% are used for assuring a real solutions (Henk Kiers
% is thanked for providing the similarity transformations)
%
% INPUTS:
% X1    : I x J matrix of data from observation one
% X2    : I x J matrix of data from observation two
% Fac   : Number of factors
%
% OUTPUTS:
% A     : Components in the row mode (I x F)
% B     : Components in the column mode (J x F)
% C     : Weights for each slab; C(1,:) are the component
%         weights for first slab such that the approximation
%         of X1 is equivalent to X1 = A*diag(C(1,:))*B.'
%

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.

% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $
% $ Version 1.03 $ Date 22. February 1999 $ Not compiled $

IsReal=0; % If complex data, complex solutions allowed.
if all(isreal(X1))&all(isreal(X2))
    IsReal=1;
end

% Find optimal bases in F x F subspace
[U,s,V]=svd(X1+X2);
U=U(:,1:F);
V=V(:,1:F);

% Reduce to an F x F dimensional subspace
S1=U'*X1*V;
S2=U'*X2*V;

% Solve eigenvalue-problem and sort according to size
[k,l]=eig(S1\S2);
l=diag(l);
ii=abs(l)>eps;
k=k(:,ii);
l=l(ii);
p=length(l);
[l,ii]=sort(l);
j=p:-1:1;
l=l(j);
l=diag(l);
k=k(:,ii(j));
k=k/norm(k);

if IsReal % Do not allow complex solutions if only reals are considered
    T1=eye(F);
    T2=eye(F);
    [rhok,argk]=complpol(k);
    [rhol,argl]=complpol(diag(l));
    j=1;
    while j<=F
        if abs(imag(l(j,j)))<.00000001  % real eigenvalue
            if abs(imag(k(1,j)))>=.00000001 % complex eigenvector
                T1(j,j)=exp(i*argk(1,j));
            end;
        end;
        if abs(imag(l(j,j)))>=.00000001  % j-th and j+1-th are complex eigenvalues
            c=argk(1,j)+argk(1,j+1);
            T1(j,j)=exp(i*c/2);
            T1(j+1,j+1)=exp(i*c/2);
            T2(j:j+1,j:j+1)=[1 1;i -i];
            j=j+1;
        end;
        j=j+1;
    end;
    k=real(k/T1/T2);
    l=T2*T1*l/T1/T2;
    l=real(diag(diag(l)));
end

C(2,:)=ones(1,F);
C(1,:)=diag(l)';
A = U*S1*k;
B=V/k';

C=(pinv(krb(B,A))*[X1(:) X2(:)]).';

function [E]=gsm(V);
%GSM orthogonalization
%
% [E]=GSM(V);
% GS   Gram-Schmidt Method for orthogonalisation
%      An orthonormal basis spanning the columns of V is returned in E.
%
%      This algorithm does not use pivoting or any other
%      stabilization scheme. For a completely safe orthogonalization
%      you should use 'ORTH()' though is may take triple the time.
%      'GSM()' is optimized for speed and requies only minimum storage
%      during iterations. No check of rank is performed on V!
%

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.


[m n]=size(V);

%Allocate space for the basis
E=zeros(m,n);

%The first basis vector is taken directly from V
s=sqrt(sum(V(:,1).^2));
E(:,1)=V(:,1)/s;

%Find the other basis vectors as orthogonals to
%the already determined basis by projection
for k=2:n,
    f=V(:,k)-E(:,1:(k-1))*(E(:,1:(k-1))'*V(:,k));
    s=sqrt(sum(f.^2));
    if s<eps,
        E(:,k)=0*f;   %set to zeros
    else
        E(:,k)=f/s;   %normalize
    end;
end;


function [X]=missmult(A,B)

%MISSMULT product of two matrices containing NaNs
%
%[X]=missmult(A,B)
%This function determines the product of two matrices containing NaNs
%by finding X according to
%     X = A*B
%If there are columns in A or B that are pur missing values,
%then there will be entries in X that are missing too.
%
%The result is standardized, that is, corrected for the lower
%number of contributing terms.
%
%Missing elements should be denoted by 'NaN's


% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.

%INBOUNDS
%REALONLY

[ia ja]=size(A);
[ib jb]=size(B);
X=zeros(ia,jb);

one_arry=ones(ia,1);
for j=1:jb,
    p=one_arry*B(:,j)';
    tmpMat=A.*p;
    X(:,j)=misssum(tmpMat')';
end;


function C=ckron(A,B)
%CKRON
% C=ckron(A,B)
%
% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.

% Claus Andersson, Jan. 1996

% Should not be compiled to overwrite ckron.mex

[mA,nA] = size(A);
[mB,nB] = size(B);

C = zeros(mA*mB,nA*nB);
if mA*nA <= mB*nB
    for i = 1:mA
        iC = 1+(i-1)*mB:i*mB;
        for j = 1:nA
            jC = 1+(j-1)*nB:j*nB;
            C(iC,jC) = A(i,j)*B;
        end
    end
else
    for i = 1:mB
        iC = i:mB:(mA-1)*mB+i;
        for j = 1:nB
            jC = j:nB:(nA-1)*nB+j;
            C(iC,jC) = B(i,j)*A;
        end
    end
end

function G=neye(Fac);
% NEYE  Produces a super-diagonal array
%
%function G=neye(Fac);
%
% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $
% $ Version 1.00 $ Date 5. Aug. 1998 $ Not compiled $
%
% This algorithm requires access to:
% 'getindxn'
%
% See also:
% 'parafac' 'maxvar3' 'maxdia3'
%
% ---------------------------------------------------------
%             Produces a super-diagonal array
% ---------------------------------------------------------
%
% G=neye(Fac);
%
% Fac      : A row-vector describing the number of factors
%            in each of the N modes. Fac must be a 1-by-N vector.
%            Ex. [3 3 3] or [2 2 2 2]

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.

N=size(Fac,2);
if N==1,
    fprintf('Specify ''Fac'' as e vector to define the order of the core, e.g.,.\n')
    fprintf('G=eyecore([2 2 2 2])\n')
end;

G=zeros(Fac(1),prod(Fac(2:N)));

for i=1:Fac(1),
    [gi,gj]=getindxn(Fac,ones(1,N)*i);
    G(gi,gj)=1;
end;

G = reshape(G,Fac);

function [b,All,MaxML]=ulsr(x,NonNeg);

%ULSR
%
% See also:
% 'unimodal' 'monreg' 'fastnnls'
%
% ------INPUT------
%
% x       is the vector to be approximated
% NonNeg  If NonNeg is one, nonnegativity is imposed
%
%
%
% ------OUTPUT-----
%
% b 	     is the best ULSR vector
% All      is containing in its i'th column the ULSRFIX solution for mode
% 	        location at the i'th element. The ULSR solution given in All
%          is found disregarding the i'th element and hence NOT optimal
% MaxML    is the optimal (leftmost) mode location (i.e. position of maximum)
%
% Reference
% Bro and Sidiropoulos, "Journal of Chemometrics", 1998, 12, 223-247.
%
%
% [b,All,MaxML]=ulsr(x,NonNeg);
% This file uses MONREG.M


% Copyright (C) 1995-2006  Rasmus Bro & Nikos Sidiropoulos
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.


% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $
%


x=x(:);
I=length(x);
xmin=min(x);
if xmin<0
    x=x-xmin;
end


% THE SUBSEQUENT
% CALCULATES BEST BY TWO MONOTONIC REGRESSIONS

% B1(1:i,i) contains the monontonic increasing regr. on x(1:i)
[b1,out,B1]=monreg(x);

% BI is the opposite of B1. Hence BI(i:I,i) holds the monotonic
% decreasing regression on x(i:I)
[bI,out,BI]=monreg(flipud(x));
BI=flipud(fliplr(BI));

% Together B1 and BI can be concatenated to give the solution to
% problem ULSR for any modloc position AS long as we do not pay
% attention to the element of x at this position


All=zeros(I,I+2);
All(1:I,3:I+2)=B1;
All(1:I,1:I)=All(1:I,1:I)+BI;
All=All(:,2:I+1);
Allmin=All;
Allmax=All;
% All(:,i) holds the ULSR solution for modloc = i, disregarding x(i),


iii=find(x>=max(All)');
b=All(:,iii(1));
b(iii(1))=x(iii(1));
Bestfit=sum((b-x).^2);
MaxML=iii(1);
for ii=2:length(iii)
    this=All(:,iii(ii));
    this(iii(ii))=x(iii(ii));
    thisfit=sum((this-x).^2);
    if thisfit<Bestfit
        b=this;
        Bestfit=thisfit;
        MaxML=iii(ii);
    end
end

if xmin<0
    b=b+xmin;
end


% Impose nonnegativity
if NonNeg==1
    if any(b<0)
        id=find(b<0);
        % Note that changing the negative values to zero does not affect the
        % solution with respect to nonnegative parameters and position of the
        % maximum.
        b(id)=zeros(size(id))+0;
    end
end

function st=stdnan(X);

%STDNAN estimate std with NaN's
%
% Estimates the standard deviation of each column of X
% when there are NaN's in X.
%
% Columns with only NaN's get a standard deviation of zero


% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.


[I,J]=size(X);

st=[];
for j=1:J
    id=find(~isnan(X(:,j)));
    if length(id)
        st=[st std(X(id,j))];
    else
        st=[st 0];
    end
end

function [Factors]=inituck(X,Fac,MthFl,IgnFl)
%INITUCK initialization of loadings
%
% function [Factors]=inituck(X,Fac,MthFl,IgnFl)
%
% This algorithm requires access to:
% 'gsm' 'fnipals' 'missmult' 'missmean'
%
% ---------------------------------------------------------
%        Initialize Factors for the Tucker3 model
% ---------------------------------------------------------
%
% [Factors]=inituck(X,Fac,MthFl,IgnFl);
% [Factors]=inituck(X,Fac);
%
% X        : The multi-way data array.
% Fac      : Vector describing the number of factors
%            in each of the N modes.
% MthFl    : Method flag indicating what kind of
%            factors you want to initiate Factors with:
%            '1' : Random values, orthogonal
%            '2' : Normalized singular vectors, orthogonal
%            '3' : SVD with successive projections
% IgnFl    : This feature is only valid with MthFl==2.
%            If specified, these mode(s) will be ignored,
%            e.g. IgnFl=[1 5] or IgnFl=[3] will
%            respectively not initialize modes one and
%            five, and mode three.
% Factors  : Contains, no matter what method, orthonormal
%            factors. This is the best general approach to
%            avoid correlated, hence ill-posed, problems.
%
% The task of this initialization program is to find acceptable
% guesses to be used as starting point in the 'TUCKER.M' program.
% Note that it IS possible to initialize the factors to have
% more columns than rows, since this may be required by some
% models. If this is required, the 'superfluos'
% columns will be random and orthogonal.
% This algorithm automatically arranges the sequence of the
% initialization to minimize time and memory consumption.
% If you get a warning from a NIPALS algorithm about convergence has
% not been reached, you can simply ignore this. With regards
% to initialization this is not important as long as the
% factors being returned are in the range of the eigensolutions.

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.

% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $
% $ Version 1.00 $ Date 24. May 1998 $ Not compiled $

format long
format compact
DimX = size(X);
X = reshape(X,DimX(1),prod(DimX(2:end)));


MissingExist=any(isnan(X(:)));

% Initialize system variables
N=size(Fac,2);
FIdx0=zeros(1,N);
FIdx1=zeros(1,N);
latest=1;
for c=1:N,
    if Fac(c)==-1,
        FIdx0(c)=0;
    else
        FIdx0(c)=latest;
        latest=latest+Fac(c)*DimX(c);
        FIdx1(c)=latest-1;
    end;
end;

% Check inputs
if ~exist('IgnFl'),
    IgnFl=[0];
end;

%Random values
if MthFl==1,
    for c=1:N,
        A=orth(rand( DimX(c) , min([Fac(c) DimX(c)]) ));
        B=[A orth(rand(DimX(c),Fac(c)-DimX(c)))];
        Factors(FIdx0(c):FIdx1(c))=B(:)';
    end;
end;

%Singular vectors
%Factors=rand(1,sum(~(Fac==-1).*DimX.*Fac)); %Matlab 4.2 compatibility
Factors=rand(1,sum((Fac~=-1).*DimX.*Fac)); %Matlab 4.2 compatibility
if MthFl==2 | MthFl==3
    
    %Remove, in a fast way, the missing values by
    %approximations as means of columns and rows
    if MissingExist,
        [i j]=find(isnan(X));
        mnx=missmean(X)/3;
        mny=missmean(X')/3;
        n=size(i,1);
        for k=1:n,
            i_=i(k);
            j_=j(k);
            X(i_,j_) = mny(i_) + mnx(j_);
        end;
        mnz=(missmean(mnx)+missmean(mny))/2;
        p=find(isnan(X));
        X(p)=mnz;
    end;
    
    [A Order]=sort(Fac);
    RedData=X;
    CurDimX=DimX;
    for k=1:N,
        c=Order(k);
        if Fac(c)>0,
            for c1=1:c-1;
                newi=CurDimX(c1+1);
                newj=prod(CurDimX)/CurDimX(c1+1);
                RedData=reshape(RedData',newi,newj);
            end;
            Op=0;
            if Op==0 & Fac(c)<=5 & (10<min(size(RedData)) & min(size(RedData))<=120),
                %Need to apply NIPALS
                A=reshape(Factors(FIdx0(c):FIdx1(c)),DimX(c),Fac(c));
                A=fnipals(RedData,min([Fac(c) DimX(c)]),A);
                B=[A orth(rand(DimX(c),Fac(c)-DimX(c)))];
                Factors(FIdx0(c):FIdx1(c))=B(:)';
                Op=1;
            end;
            if Op==0 & (120<min(size(RedData)) & min(size(RedData))<Inf),
                %Need to apply Gram-Schmidt
                C=RedData*RedData';
                A=reshape(Factors(FIdx0(c):FIdx1(c)),DimX(c),Fac(c));
                for i=1:3,
                    A=gsm(C*A);
                end;
                B=[A orth(rand(DimX(c),Fac(c)-DimX(c)))];
                Factors(FIdx0(c):FIdx1(c))=B(:)';
                Op=1;
            end;
            if Op==0 & (0<min(size(RedData)) & min(size(RedData))<=120),
                %Small enough to apply SVD
                [U S A]=svd(RedData',0);
                A=A(:,1:min([Fac(c) DimX(c)]));
                B=[A orth(rand(DimX(c),Fac(c)-DimX(c)))];
                Factors(FIdx0(c):FIdx1(c))=B(:)';
                Op=1;
            end;
            CurDimX(c)=min([Fac(c) DimX(c)]);
            RedData=A'*RedData;
            %Examine if re-ordering is necessary
            if c~=1,
                for c1=c:N,
                    if c1~=N,
                        newi=CurDimX(c1+1);
                        newj=prod(CurDimX)/newi;
                    else
                        newi=CurDimX(1);
                        newj=prod(CurDimX)/newi;
                    end;
                    RedData=reshape(RedData',newi,newj);
                end;
            end;
        end;
    end;
end;
format

% Convert to new format
clear ff
id1 = 0;

for i = 1:length(DimX)
    
    if Fac(i)~=-1
        id2 = sum(DimX(1:i).*Fac(1:i).*(Fac(1:i)~=-1));
        ff{i} = reshape(Factors(id1+1:id2),DimX(i),Fac(i));
        id1 = id2;
    else
        ff{i}=[];
    end
end
Factors = ff;

function [mm]=misssum(X,def)
%MISSSUM sum of a matrix X with NaN's
%
%[mm]=misssum(X,def)
%
%This function calculates the sum of a matrix X.
%X may hold missing elements denoted by NaN's which
%are ignored.
%
%The result is standardized, that is, corrected for the lower
%number of contributing terms.
%
%Check that for no column of X, all values are missing

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.

%Insert zeros for missing, correct afterwards
missidx = isnan(X);
i = find(missidx);
if ~isempty(i),
    X(i) = zeros(size(i));
end;

%Find the number of real(non-missing objects)
if min(size(X))==1,
    n_real=length(X)-sum(missidx);
    weight=length(X);
else
    n_real=size(X,1)-sum(missidx);
    weight=size(X,1);
end

i=find(n_real==0);
if isempty(i) %All values are real and can be corrected
    mm=weight*sum(X)./n_real;
else %There are columns with all missing, insert missing
    n_real(i)=1;
    mm=weight*sum(X)./n_real;
    mm(i)=i + NaN;
end

function mm=missmean(X)
%MISSMEAN Mean of a matrix X with NaN's.
%
%I/O: [mm]=missmean(X)
%

% Copyright Eigenvector Research, Inc. 2001-2011
%Licensee shall not re-compile, translate or convert "M-files" contained
% in PLS_Toolbox for use with any software other than MATLAB®, without
% written permission from Eigenvector Research, Inc.

%Insert zeros for missing, correct afterwards
missidx = isnan(X);
i = find(missidx);
X(i) = 0;

%Find the number of real(non-missing objects)
if min(size(X))==1,
    n_real=length(X)-sum(missidx);
else
    n_real=size(X,1)-sum(missidx);
end

i=find(n_real==0);
if isempty(i) %All values are real and can be corrected
    mm=sum(X)./n_real;
else %There are columns with all missing, insert missing
    n_real(i)=1;
    mm=sum(X)./n_real;
    mm(i)=i + NaN;
end


function [G]=calcore(X,Factors,Options,O,MissingExist);

%CALCORE Calculate the Tucker core
%
%
% [G]=calcore(X,Factors,Options);
% [G]=calcore(X,Factors);
%
% This algorithm applies to the general N-way case, so
% the unfolded X can have any number of dimensions. The principles of
% 'projections' and 'systematic unfolding methodology (SUM)' are used
% in this algorithm so orthogonality is required.
% This algorithm can handle missing values in X and
% also allows for TUCKER2 models using the an empty matrix in the
% corresponding cell of Factors.
% The variable 'Factors' must contain the stringed-out factors.

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.

format compact
format long

DimX = size(X);
X = reshape(X,DimX(1),prod(DimX(2:end)));
ff = [];
for f=1:length(Factors)
    ff=[ff;Factors{f}(:)];
    Fac(f)=size(Factors{f},2);
    if isempty(Factors{f}) % 'Tucker2' - i.e. no compression in that mode
        Fac(f) = -1;
    end
end
Factors = ff;

% Initialize system variables
if length(Fac)==1,
    Fac=Fac*ones(size(DimX));
end;

Fac_orig=Fac;
i=find(Fac==-1);
Fac(i)=zeros(1,length(i));
N=size(Fac,2);
FIdx0=zeros(1,N);
FIdx1=zeros(1,N);
if ~exist('MissingExist')
    if sum(isnan(X(:)))>0,
        MissingExist=1;
    else
        MissingExist=0;
    end;
end;
FIdx0=cumsum([1 DimX(1:N-1).*Fac(1:N-1)]);
FIdx1=cumsum([DimX.*Fac]);
if ~exist('O') | isempty(O),
    O=1;
end;


if O, %means orthogonality
    CurDimX=DimX;
    RedData=X;
    for c=1:N,
        
        if Fac_orig(c)==-1,
            kthFactor=eye(DimX(c));
            CurDimX(c)=DimX(c);
        else
            kthFactor=reshape(Factors(FIdx0(c):FIdx1(c)),DimX(c),Fac(c));
            CurDimX(c)=Fac(c);
        end;
        if MissingExist
            RedData=missmult(kthFactor',RedData);
        else
            RedData=kthFactor'*RedData;
        end;
        
        if c~=N,
            newi=CurDimX(c+1);
            newj=prod(CurDimX)/CurDimX(c+1);
        else
            newi=CurDimX(1);
            newj=prod(CurDimX)/CurDimX(1);
        end;
        
        RedData=reshape(RedData',newi,newj);
    end;
    G=RedData;
else %oblique factors
    
    LMatTmp=1;
    if Fac_orig(1)==-1,
        LMatTmp=eye(DimX(c));
    else
        LMatTmp=reshape(Factors(FIdx0(1):FIdx1(1)),DimX(1),Fac(1));
    end;
    
    RMatTmp=1;
    for c=2:N,
        if Fac_orig(c)==-1,
            kthFactor=eye(DimX(c));
        else
            kthFactor=reshape(Factors(FIdx0(c):FIdx1(c)),DimX(c),Fac(c));
        end;
        RMatTmp=ckron(kthFactor',RMatTmp);
    end;
    
    if MissingExist
        RedData=missmult(pinv(LMatTmp),X);
        RedData=missmult(RedData,pinv(RMatTmp));
    else
        RedData=LMatTmp\X;
        RedData=RedData/RMatTmp;
    end;
    
    G=RedData;
    
end;

for i = 1:length(Fac)
    if Fac(i)==0
        Fac(i) = DimX(i);
    end
end
G = reshape(G,Fac);

return



function [F,SS]=nonneg(X,M,F,Options);
%NONNEG1 alternative to NNLS
%
%function [F,SS]=nonneg(X,M,F,Options);
%
% 'nonneg.m'
% $ Version 0.01 $ Date 6. Aug. 1997 $ Not compiled $
%
% This algorithm requires access to:
% ''
%
% ---------------------------------------------------------
%             Fast! Non-negativity Regression
% ---------------------------------------------------------
%
% [F,SS]=nonneg(X,M,F,Options);
% [F,SS]=nonneg(X,M);
%
% X        : Matrix of regressors.
% M        : Matrix of regressands.
% F        : The non-negativity constrained solution matrix.
% Options  : Type 'help Options' at the prompt.
%
% Given X and M this algorithm solves for the optimal
% F in a least squares sense, using that
%      X = F*M
% in the problem
%      min ||X-F*M||, s.t. F>=0, for given X and M.
%
% This version does not accept missing values.

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.

format long
format compact
Show=0;

SS=[];

[aX bX]=size(X);
[aM bM]=size(M);
aF=aX;
bF=aM;
w=bF;

Xc=X;
W=eye(w);
itimax=100;
itomax=100;
SSimax=1e-10;
SSomax=1e-12;
SSiOld=realmax;
SSoOld=realmax;
%Initialize F
if ~exist('F')
    F=X*M'/(M*M');
    FOld=F;
    I=find(F<0);
    FOld(I)=0;
    SSiOld=sum(sum( (F-FOld).^2 ));
    F=FOld;
end;
FOld=F;
SSX=sum(sum(X.^2));
MMT=M*M';
XMT=X*M';
InvMMT=1./diag(MMT);

ito=0;
convo=0;
while ~convo,
    ito=ito+1;
    convi=0;
    iti=0;
    while ~convi,
        iti=iti+1;
        
        %Iterate on variables for non-negativity
        for i=1:w,
            W(i,i)=0;
            f=XMT(:,i)-F*W*MMT(:,i);
            f=InvMMT(i)*f;
            I=find(f<0);
            if ~isempty(I),
                %f(I)=0;
                f(I)=zeros(1,length(I));
            end;
            F(:,i)=f;
            W(i,i)=1;
        end;
        
        %Estimate error now
        SSi=sum(sum( (F-FOld).^2 ));
        if SSi < sum(sum( FOld.^2 ))*SSimax | iti>itimax,
            convi=1;
        end;
        FOld=F;
        
    end;
    
    %Estimate error on the transformed LS problem
    SSo=sum(sum( (XMT-F*MMT).^2 ));
    if (SSoOld-SSo)/SSX<SSomax,
        convo=1;
    end;
    if ito>itomax,
        convo=1;
    end;
    SSoOld=SSo;
end;
format


function [Y,D]=normit(X)
%NORMIT normalize
%
%[Y,D]=normit(X);
%Normalizes X, such that X  = Y *D
%It also applies that    X' = D'*Y'
%Always give X column wise


% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.

[a b]=size(X);

Y=zeros(a,b);

SS=sqrt(sum(X.^2));
for i=1:b,
    Y(:,i)=X(:,i)./SS(i);
end;

if nargout==2,
    D=(Y'*Y)\(Y'*X);
end;

function B=unimodal(X,Y,Bold)

%UNIMODAL unimodal regression
%
% Solves the problem min|Y-XB'| subject to the columns of
% B are unimodal and nonnegative. The algorithm is iterative
% If an estimate of B (Bold) is given only one iteration is given, hence
% the solution is only improving not least squares
% If Bold is not given the least squares solution is estimated
%
% I/O B=unimodal(X,Y,Bold)
%
% Reference
% Bro and Sidiropoulos, "Journal of Chemometrics", 1998, 12, 223-247.



% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.


if nargin==3
    B=Bold;
    F=size(B,2);
    for f=1:F
        y=Y-X(:,[1:f-1 f+1:F])*B(:,[1:f-1 f+1:F])';
        beta=pinv(X(:,f))*y;
        B(:,f)=ulsr(beta',1);
    end
else
    F=size(X,2);
    maxit=100;
    B=randn(size(Y,2),F);
    Bold=2*B;
    it=0;
    while norm(Bold-B)/norm(B)>1e-5&it<maxit
        Bold=B;
        it=it+1;
        for f=1:F
            y=Y-X(:,[1:f-1 f+1:F])*B(:,[1:f-1 f+1:F])';
            beta=pinv(X(:,f))*y;
            B(:,f)=ulsr(beta',1);
        end
    end
    if it==maxit
        disp([' UNIMODAL did not converge in ',num2str(maxit),' iterations']);
    end
end

function [v]=nsetdiff(A,B);
%NSETDIFF
%
%[v]=nsetdiff(A,B);
%Slow setdiff by CA, 1998

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.

AS=sort(A);
len_AS=length(AS);
for i=1:len_AS-1,
    for j=i+1:len_AS,
        if AS(i)==AS(j),
            AS(i)=NaN;
        end;
    end
end;
I=find(isnan(AS));
if ~isempty(I)
    AS(I)=[];
end;


BS=sort(B);
len_BS=length(BS);
for i=1:len_BS-1,
    for j=i+1:len_BS,
        if BS(i)==BS(j),
            BS(i)=NaN;
        end;
    end
end;
I=find(isnan(BS));
if ~isempty(I)
    BS(I)=[];
end;


len_AS=length(AS);
len_BS=length(BS);
if len_AS >= len_BS
    for i=1:len_AS,
        for j=1:len_BS,
            if AS(i)==BS(j),
                AS(i)=NaN;
            end;
        end;
    end;
    I=find(isnan(AS));
    if ~isempty(I)
        AS(I)=[];
    end;
    v=AS;
else
    for i=1:len_BS,
        for j=1:len_AS,
            if BS(i)==AS(j),
                BS(i)=NaN;
            end;
        end;
    end;
    I=find(isnan(BS));
    if ~isempty(I)
        BS(I)=[];
    end;
    v=BS;
end;



function [R,Y]=complpol(X);
%COMPLPOL
% produces radius and angle for matrix of complex data
R=(real(X).^2+imag(X).^2).^.5;
A=X./R;
Y=real(log(A)/i);

function [i,j]=getindxn(R,Idx);
%GETINDXN
%
%[i,j]=getindxn(R,Idx)


% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.

l=size(Idx,2);

i=Idx(1);
j=Idx(2);

if l==3,
    j = j + R(2)*(Idx(3)-1);
else
    for q = 3:l,
        j = j + prod(R(2:(q-1)))*(Idx(q)-1);
    end;
end;

function [b,B,AllBs]=monreg(x);

%MONREG monotone regression
%
% See also:
% 'unimodal' 'monreg' 'fastnnls'
%
%
% MONTONE REGRESSION
% according to J. B. Kruskal 64
%
% b     = min|x-b| subject to monotonic increase
% B     = b, but condensed
% AllBs = All monotonic regressions, i.e. AllBs(1:i,i) is the
%         monotonic regression of x(1:i)
%
% Reference
% Bro and Sidiropoulos, "Journal of Chemometrics", 1998, 12, 223-247.
%
% [b,B,AllBs]=monreg(x);

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.

I=length(x);
if size(x,2)==2
    B=x;
else
    B=[x(:) ones(I,1)];
end

AllBs=zeros(I,I);
AllBs(1,1)=x(1);
i=1;
while i<size(B,1)
    if B(i,1)>B(min(I,i+1),1)
        summ=B(i,2)+B(i+1,2);
        B=[B(1:i-1,:);[(B(i,1)*B(i,2)+B(i+1,1)*B(i+1,2))/(summ) summ];B(i+2:size(B,1),:)];
        OK=1;
        while OK
            if B(i,1)<B(max(1,i-1),1)
                summ=B(i,2)+B(i-1,2);
                B=[B(1:i-2,:);[(B(i,1)*B(i,2)+B(i-1,1)*B(i-1,2))/(summ) summ];B(i+1:size(B,1),:)];
                i=max(1,i-1);
            else
                OK=0;
            end
        end
        bInterim=[];
        for i2=1:i
            bInterim=[bInterim;zeros(B(i2,2),1)+B(i2,1)];
        end
        No=sum(B(1:i,2));
        AllBs(1:No,No)=bInterim;
    else
        i=i+1;
        bInterim=[];
        for i2=1:i
            bInterim=[bInterim;zeros(B(i2,2),1)+B(i2,1)];
        end
        No=sum(B(1:i,2));
        AllBs(1:No,No)=bInterim;
    end
end

b=[];
for i=1:size(B,1)
    b=[b;zeros(B(i,2),1)+B(i,1)];
end


function lambda_b=maxlam(X)

[I,J]=size(X);
% Ibar is the expected (average) co-cluster size (membership) in the row
% mode. E.g., if I=50 but we know that about 10 will belong to one
% co-cluster, we use Ibar = 10.
Ibar = max(10,round(I/10));
acn=zeros(1,J);
for j=1:J,
    acn(j)=norm(X(:,j),2);
end
mcn=max(acn);
lambda_bstar = 2*max(max(abs(X)))*Ibar*mcn;
lambda_b = 0.05*lambda_bstar;

function [lamA,lamB,go_on,direction]=adjustlambda(lamA,lamB,A,B,lambda_Amax,lambda_Bmax,scale);

direction=[0 0];
showit=1;
go_on=0;
lamAold=lamA;lamBold=lamB;
if any(sum(A==0)==size(A,1)) % If any column has all zeros
    % Then lower a bit (on a log scale)
    lamA = exp(log(1+lamA)*.8^scale)-1;
    if showit;disp('Lowering lambda on A');end
    go_on=1;
    direction(1)=-1; % Decreasing
elseif max(sum(A==0))<(size(A,1)/.8) % If no column has many zeros
    if ~any(sum(A==0)==size(A,1))
    % Then increase a bit (on a log scale)
    lamA = exp(log(1+lamA)*1.2^scale)-1;
    if showit;disp('Increasing lambda on A');end
    go_on=1;
    direction(1)=1; % Increasing
    end
end

if any(sum(B==0)==size(B,1)) % If any column has all zeros
    % Then lower a bit (on a log scale)
    lamB = exp(log(1+lamB)*.8^scale)-1;
    if showit;disp('Lowering lambda on B');end
    go_on=1;
    direction(2)=-1; % Decreasing
 elseif any(sum(B==0))<(size(B,1)/.8) % If no column has many zeros
    if ~any(sum(B==0)==size(B,1))
    % Then increase a bit (on a log scale)
    lamB = exp(log(1+lamB)*1.2^scale)-1;
    if showit;disp('Increasing lambda on B');end
    go_on=1;
    direction(2)=1; % Increasing
    end
end
if isnan(lamA)|isnan(lamB);
    save this
    error('SSSHSH');
end
    