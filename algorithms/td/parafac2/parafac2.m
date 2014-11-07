function [A,H,C,P,fit,AddiOutput]=parafac2(X,F,Constraints,Options,A,H,C,P);

% $ Version 1.01 $ Date 28. December 1998 $ Not compiled $ RB

% $ Version 1.02 $ Date 31. March    1999 $ Added X-validation and added function $ Not compiled $ RB
% $ Version 1.03 $ Date 20. April    1999 $ Cosmetic changes $ Not compiled $ RB
% $ Version 1.04 $ Date 25. April    1999 $ Cosmetic changes $ Not compiled $ RB
% $ Version 1.05 $ Date 18. May      1999 $ Added orthogonality constraints $ Not compiled $ RB
% $ Version 1.06 $ Date 14. September1999 $ Changed helpfile $ Not compiled $ RB
% $ Version 1.07 $ Date 20. October  1999 $ Added unimodality $ Not compiled $ RB
% $ Version 1.08 $ Date 27. March    2000 $ Optimized handling of missing dat $ Not compiled $ RB
%

% This M-file and the code in it belongs to the holder of the 
% copyrights and is made public under the following constraints:
% It must not be changed or modified and code cannot be added.
% The file must be regarded as read-only. Furthermore, the
% code can not be made part of any toolbox or similar.
% In case of doubt, contact the holder of the copyrights.

%
% Rasmus Bro
% Chemometrics Group, Food Technology
% Department of Food and Dairy Science
% Royal Veterinary and Agricultutal University
% Rolighedsvej 30, DK-1958 Frederiksberg, Denmark
% Phone  +45 35283296
% Fax    +45 35283245
% E-mail rb@kvl.dk
%

%     ___________________________________________________
%
%                  THE PARAFAC2 MODEL
%     ___________________________________________________
% 
%

% Algorithm to fit the PARAFAC2 model which is an advanced variant of the 
% normal PARAFAC1 model. It handles slab-wise deviations between components
% in one mode as long as the cross-product of the components stays 
% reasonably fixed. This can be utilized for modeling chromatographic 
% data with retention time shifts, modeling certain batch data of 
% varying length etc. See Bro, Kiers & Andersson, Journal of Chemometrics,
% 1999, 13, 295-309 for details on application and Kiers, ten Berge & 
% Bro, Journal of Chemometrics, 1999, 13, 275-294, for details on the algorithm
% 
%
% The PARAFAC2 model is given
% 
% Xk = A*Dk*(Pk*H)' + Ek, k = 1, .., K
% 
% Xk is a slab of data (I x J) in which J may actually vary with K. K 
% is the number of slabs. A (I x F) are the scores or first-mode loadings. Dk 
% is a diagonal matrix that holds the k'th row of C in its diagonal. C 
% (K x F) is the third mode loadings, H is an F x F matrix, and Pk is a
% J x F orthogonal matrix (J may actually vary from k to k. The output here
% is given as a cell array of size J x F x K. Thus, to get e.g. the second P
% write P(:,:,2), and to get the estimate of the second mode loadings at this
% second frontal slab (k = 2), write P(:,:,2)*H. The matrix Ek holds the residuals.
% 
% INPUT
% 
% X
%   Holds the data.
%   If all slabs have similar size, X is an array:
%      X(:,:,1) = X1; X(:,:,2) = X2; etc.  
%   If the slabs have different size X is a cell array (type <<help cell>>)
%      X{1} = X1; X{2} = X2; etc.
%   If you have your data in an 'unfolded' two-way array of size
%   I x JK (the three-way array is I x J x K), then simply type
%   X = reshape(X,[I J K]); to convert it to an array.
%
% F
%   The number of components to extract
% 
% Constraints
%   Vector of length 2. The first element defines constraints
%   imposed in the first mode, the second defines contraints in
%   third mode (the second mode is not included because constraints
%   are not easily imposed in this mode)
% 
%   If Constraints = [a b], the following holds. If 
%   a = 0 => no constraints in the first mode
%   a = 1 => nonnegativity in the first mode
%   a = 2 => orthogonality in the first mode
%   a = 3 => unimodality (and nonnegativity) in the first mode
%   same holds for b for the third mode
%
% Options
%   An optional vector of length 3
%   Options(1) Convergence criterion
%            1e-7 if not given or given as zero
%   Options(2) Maximal iterations
%            default 2000 if not given or given as zero
%   Options(3) Initialization method
%            A rather slow initialization method is used per default
%            but it pays to investigate in avoiding local minima.
%            Experience may point to faster methods (set Options(3)
%            to 1 or 2). You can also change the number of refits etc.
%            in the beginning of the m-file
%            0 => best of 10 runs of maximally 80 iterations (default)
%            1 => based on SVD
%            2 => random numbers
%   Options(4) Cross-validation
%            0 => no cross-validation
%            1 => cross-validation splitting in 7 segments
%   Options(5) show output
%            0 => show standard output on screen
%            1 => hide all output to screen
%
% AUXILIARY
% - Missing elements: Use NaN for missing elements
% - You can input initial values by using the input argument
%           (X,F,Constraints,Options,A,H,C,P);
%
% OUTPUT
% See right above INPUT
% 
% I/O
% 
% Demo
% parafac2('demo')
% 
% Short 
% [A,H,C,P]=parafac2(X,F);
%
% Long
% [A,H,C,P,fit]=parafac2(X,F,Constraints,Options);
%
% Copyright
% Rasmus Bro
% KVL, DK, 1998
% rb@kvl.dk
%
% Reference to algorithm
% Bro, Kiers & Andersson, PARAFAC2 - Part II. Modeling chromatographic 
% data with retention time shifts, Journal of Chemometrics, 1999, 13, 295-309

% TO DO:
% Set the algorithm to handle fixed modes as in PARALIN
% Make it N-way
% Incorporate ulsr


if nargin==0

   disp(' ')

   disp(' ')

   disp(' THE PARAFAC2 MODEL')

   disp(' ')

   disp(' Type <<help parafac2>> for more info')

   disp('  ')
   disp(' I/O ')
   disp(' [A,H,C,P]=parafac2(X,F);')
   disp(' ')

   disp(' Or optionally')
   disp(' ')

   disp(' [A,H,C,P,fit]=parafac2(X,F,Constraints,Options);')

   disp(' ')

   disp(' Options=[Crit MaxIt Init Xval Show]')

   disp(' ')

   disp(' ')

   return;

 elseif nargin<2&~all(X=='demo')

    error(' The inputs X and F must be given')

 end

 
  if isstr(X) & all(X=='demo')
    F=3;
    n=1:30;
    disp(' ')
    disp(' %%%%% PARAFAC2 DEMO %%%%%%')
    disp(' ')
    disp(' Generating simulated data')
    disp(' Note that the second mode loadings change from slab to slab')
    disp(' hence the ordinary PARAFAC model is not valid')
    disp(' ')
    subplot(2,2,1)
    A=[exp(-((n-15)/5).^2);exp(-((n-1)/10).^2);exp(-((n-21)/7).^2)]';
    plot(A),title(' First mode loadings')
    subplot(2,2,2)
    C=rand(4,3);
    plot(C),title(' Third mode loadings')
    H=orth(orth(rand(F))');
    P=[];X=[];
    for i=1:size(C,1),
       subplot(2,4,4+i)
       P(:,:,i)=orth(rand(7,F));
       plot(P(:,:,i)*H),eval(['title([''2. mode, k = '',num2str(i)])'])
    end,
    disp(' Press key to continue'),pause
    for i=1:size(C,1),
       X(:,:,i)=A*diag(C(i,:))*(P(:,:,i)*H)';
    end,
    
    X = X + randn(size(X))*.01;
    disp(' Adding one percent noise and fitting model')
    disp(' Several initial models will be fitted and the best used')
    [a,h,c,p]=parafac2(X,F);
    
    disp(' ')
    disp(' Results shown in plot')
    subplot(2,2,1)
    plot(A*diag(sum(A).^(-1)),'r'),
    hold on,
    plot(a*diag(sum(a).^(-1)),'g'),title(' First mode (red true,green estimated)')
    hold off
    subplot(2,2,2)
    plot(C*diag(sum(C).^(-1)),'r')
    hold on,
    plot(c*diag(sum(c).^(-1)),'g'),title(' Third mode (red true,green estimated)')
    hold off
    for i=1:size(C,1),
       subplot(2,4,4+i)
       ph=P(:,:,i)*H;
       plot(ph*diag(sum(ph).^(-1)),'r'),
       hold on
       ph=p{i}*h;
       plot(ph*diag(sum(ph).^(-1)),'g'),
       eval(['title([''2. mode, k = '',num2str(i)])'])
       hold off
    end,
    return;
   
end


ShowFit  = 1000; % Show fit every 'ShowFit' iteration
NumRep   = 10; %Number of repetead initial analyses
NumItInRep = 80; % Number of iterations in each initial fit
if ~(length(size(X))==3|iscell(X))
   error(' X must be a three-way array or a cell array')
end
%set random number generators
randn('state',sum(100*clock));
rand('state',sum(100*clock));

if nargin < 4
  Options = zeros(1,5);
end
if length(Options)<5
   Options = Options(:);
   Options = [Options;zeros(5-length(Options),1)];
end

% Convergence criterion
if Options(1)==0
   ConvCrit = 1e-7;
else
   ConvCrit = Options(1);
end
if Options(5)==0
   disp(' ')
   disp(' ')
   disp([' Convergence criterion        : ',num2str(ConvCrit)])
end

% Maximal number of iterations 
if Options(2)==0
   MaxIt = 2000;
else
   MaxIt = Options(2);
end

% Initialization method
initi = Options(3);

if nargin<3
  Constraints = [0 0];
end
if length(Constraints)~=2
   Constraints = [0 0];
   disp(' Length of Constraints must be two. It has been set to zeros')
end
% Modify to handle GPA (Constraints = [10 10]);
if Constraints(2)==10
   Constraints(1)=0;
   ConstB = 10;
else
   ConstB = 0;
end


ConstraintOptions=[ ...
   'Fixed                     ';...
   'Unconstrained             ';...
   'Non-negativity constrained';...
   'Orthogonality constrained ';...
   'Unimodality constrained   ';...
   'Not defined               ';...
   'Not defined               ';...
   'Not defined               ';...
   'Not defined               ';...
   'Not defined               ';...
   'Not defined               ';...
   'GPA                       '];
   

if Options(5)==0
   disp([' Maximal number of iterations : ',num2str(MaxIt)])
   disp([' Number of factors            : ',num2str(F)])
   disp([' Loading 1. mode, A           : ',ConstraintOptions(Constraints(1)+2,:)])
   disp([' Loading 3. mode, C           : ',ConstraintOptions(Constraints(2)+2,:)])
   disp(' ')
end


% Make X a cell array if it isn't
if ~iscell(X)
  for k = 1:size(X,3)
    x{k} = X(:,:,k);
  end
  X = x;
  clear x
end
I = size(X{1},1);
K = max(size(X));

% CROSS-VALIDATION
if Options(4)==1
   Opt = Options;
   Opt(4) = 0;
   splits = 7;
   while rem(I,splits)==0 % Change the number of segments if 7 is a divisor in prod(size(X))
      splits = splits + 2;
   end
   AddiOutput.NumberOfSegments = splits;
   if Options(5)==0
      disp(' ')
      disp([' Cross-validation will be performed using ',num2str(splits),' segments'])
      disp([' and using from 1 to ',num2str(F),' components'])
      XvalModel = [];
   end
   SS = [];
   for f = 1:F
      Arep = [];Hrep = [];Crep = [];clear Prep;
      for s = 1:splits
         Xmiss = X;
         for k = 1:K 
            Xmiss{k}(s:splits:end)=NaN;
         end
         [a,h,c,p]=parafac2(Xmiss,f,Constraints,Opt);
         Arep(:,:,s)=a;Hrep(:,:,s)=h;Crep(:,:,s)=c;Prep(s,:)=p;
         M=[];
         for k = 1:K
            m    = a*diag(c(k,:))*(p{k}*h)';
            M{k} = m;
         end
         XvalModel{f} = M;
      end
      AddiOutput.XvalModels=XvalModel;
      ss = 0;
      for k = 1:K 
         x = X{k};m = M{k};
         ss = ss + sum((x(s:splits:end)-m(s:splits:end)).^2);
      end
      SS = [SS ss];
      AddiOutput.SS = SS;
      AddiOutput.A_xval{f}=Arep;
      AddiOutput.H_xval{f}=Hrep;
      AddiOutput.C_xval{f}=Crep;
      AddiOutput.P_xval{f}=Prep;
   end

      clf
      plot([1:F],SS),title(' Residual sum-squares - cross-validation')
      xlabel('Number of components')
      disp(' ')
      disp(' The total model has NOT been fitted.')
      disp(' You must refit the model with the number of ')
      disp(' components you judge necessary.')
      disp(' ')
      disp(' You can also check the outputted struct array')
      disp(' It contains loadings estimated from different')
      disp(' subsets and stability of subsets indicates validity.')
      disp(' (e.g. if name of struct array is Output then the file')
      disp(' AX=Output.A_xval{3}is a three-way array holding all A')
      disp(' loadings estimated with 3 components. AX(:,:,1) is the ')
      disp(' estimate of A obtained from the first subset etc.')
      
      [a,b]=min(SS);
      figure
      subplot(2,1,1)
      a=AddiOutput.A_xval{b};
      for i=2:splits
         plot(a(:,:,i),'r'),hold on
      end
      title([' A resampled during Xval for ',num2str(b),' comp.'])
      hold off
      
      subplot(2,1,2)
      c = AddiOutput.C_xval{b};
      for i=2:splits
         plot(c(:,:,i),'r'),hold on
      end
      title([' C resampled during Xval for ',num2str(b),' comp.'])
      hold off
      return;
end

% Find missing and replace with average 
MissingElements = 0;
MissNum=0;AllNum=0;
for k = 1:K
   x=X{k};
   miss = sparse(isnan(x));
   MissingOnes{k} = miss;
   if any(miss(:))
     MissingElements = 1;
     % Replace missing with mean over slab (not optimal but what the heck)
     % Iteratively they'll be replaced with model estimates
     x(find(miss)) = mean(x(find(~miss)));
     X{k} = x;
     MissNum = MissNum + prod(size(find(miss)));
     AllNum = AllNum + prod(size(x));
   end
end
if MissingElements
   if Options(5)==0
      PercMiss = 100*MissNum/AllNum;
      RoundedOf = .1*round(PercMiss*10);
      disp([' Missing data handled by EM   : ',num2str(RoundedOf),'%'])
   end
end
clear x

% Initialize by ten small runs
if nargin<5
   if initi==0
      if Options(5)==0
         disp([' Use best of ',num2str(NumRep)]) 
         disp(' initially fitted models')
      end
      Opt = Options;
      Opt = Options(1)/20;
      Opt(2) = NumItInRep; % Max NumItInRep iterations
      Opt(3) = 1;  % Init with SVD
      Opt(4) = 0;
      Opt(5) = 1;
      [A,H,C,P,bestfit]=parafac2(X,F,Constraints,Opt);
      AllFit = bestfit;
      for i = 2:NumRep
         Opt(3) = 2;   % Init with random
         [a,h,c,p,fit]=parafac2(X,F,Constraints,Opt);
         AllFit = [AllFit fit];
         if fit<bestfit
            A=a;H=h;C=c;P=p;
            bestfit = fit;
         end
      end
      AddiOutput.AllFit = AllFit;
      if Options(5)==0
         for ii=1:length(AllFit)
            disp([' Initial Model Fit            : ',num2str(AllFit(ii))])
         end
      end
      % Initialize by SVD
   elseif initi==1
      if Options(5)==0
         disp(' SVD based initialization')
      end
      XtX=X{1}*X{1}';
      for k = 2:K
         XtX = XtX + X{k}*X{k}';
      end
      [A,s,v]=svd(XtX,0);  
      A=A(:,1:F);
      C=ones(K,F)+randn(K,F)/10;
      H = eye(F);
   elseif initi==2
      if Options(5)==0
         disp(' Random initialization')
      end
      A = rand(I,F);
      C = rand(K,F);
      H = eye(F);
   else
      error(' Options(2) wrongly specified')
   end
end

if initi~=1
   XtX=X{1}*X{1}'; % Calculate for evaluating fit (but if initi = 1 it has been calculated)
   for k = 2:K
      XtX = XtX + X{k}*X{k}';
   end
end  
fit    = sum(diag(XtX));
oldfit = fit*2;
fit0   = fit;
it     = 0;
Delta = 1;

if Options(5)==0
   disp(' ')
   disp(' Fitting model ...')
   disp(' Loss-value      Iteration     %VariationExpl')
end

% Iterative part
while abs(fit-oldfit)>oldfit*ConvCrit & it<MaxIt & fit>1000*eps
    oldfit = fit;
    it   = it + 1;
    
    % Update P
    for k = 1:K
      Qk       = X{k}'*(A*diag(C(k,:))*H');
      P{k}     = Qk*psqrt(Qk'*Qk);
    %  [u,s,v]  = svd(Qk.');P{k}  = v(:,1:F)*u(:,1:F)';
      Y(:,:,k) = X{k}*P{k};
    end
    
    % Update A,H,C using PARAFAC-ALS
    [A,H,C,ff]=parafac(reshape(Y,I,F*K),[I F K],F,1e-4,[Constraints(1) ConstB Constraints(2)],A,H,C,5);
    [fit,X] = pf2fit(X,A,H,C,P,K,MissingElements,MissingOnes);
      
    % Print interim result
    if rem(it,ShowFit)==0|it == 1
       if Options(5)==0
          fprintf(' %12.10f       %g        %3.4f \n',fit,it,100*(1-fit/fit0));
          subplot(2,2,1)
          plot(A),title('First mode')
          subplot(2,2,2)
          plot(C),title('Third mode')
          subplot(2,2,3)
          plot(P{1}*H),title('Second mode (only first k-slab shown)')
          drawnow
       end
    end

end

if rem(it,ShowFit)~=0 %Show final fit if not just shown
   if Options(5)==0
      fprintf(' %12.10f       %g        %3.4f \n',fit,it,100*(1-fit/fit0));
   end
end



function [fit,X]=pf2fit(X,A,H,C,P,K,MissingElements,MissingOnes);

   % Calculate fit and impute missing elements from model

   fit = 0;
   for k = 1:K
     M   = A*diag(C(k,:))*(P{k}*H)';
     % if missing values replace missing elements with model estimates
     if nargout == 2 
       if any(MissingOnes{k})
         x=X{k};
         x(find(MissingOnes{k})) = M(find(MissingOnes{k}));
         X{k} = x;
       end
     end
     fit = fit + sum(sum(abs (X{k} - M ).^2));
   end


function X = psqrt(A,tol)

   % Produces A^(-.5) even if rank-problems

   [U,S,V] = svd(A,0);
   if min(size(S)) == 1
     S = S(1);
   else
     S = diag(S);
   end
   if (nargin == 1)
     tol = max(size(A)) * S(1) * eps;
   end
   r = sum(S > tol);
   if (r == 0)
     X = zeros(size(A'));
   else
     S = diag(ones(r,1)./sqrt(S(1:r)));
     X = V(:,1:r)*S*U(:,1:r)';
  end
  
  
function [A,B,C,fit,it] = parafac(X,DimX,Fac,crit,Constraints,A,B,C,maxit,DoLineSearch);

% Complex PARAFAC-ALS
% Fits the PARAFAC model Xk = A*Dk*B.' + E
% where Dk is a diagonal matrix holding the k'th
% row of C.
%
% Uses on-the-fly projection-compression to speed up 
% the computations. This requires that the first mode 
% is the largest to be effective
% 
% INPUT
% X          : Data
% DimX       : Dimension of X
% Fac        : Number of factors
% OPTIONAL INPUT
% crit       : Convergence criterion (default 1e-6)
% Constraints: [a b c], if e.g. a=0 => A unconstrained, a=1 => A nonnegative
% A,B,C      : Initial parameter values
%
% I/O
% [A,B,C,fit,it]=parafac(X,DimX,Fac,crit,A,B,C);
%
% Copyright 1998
% Rasmus Bro
% KVL, Denmark, rb@kvl.dk

% Initialization
if nargin<9
  maxit   = 2500;      % Maximal number of iterations
end
showfit = pi;         % Show fit every 'showfit'th iteration (set to pi to avoid)

if nargin<4
  crit=1e-6;
end

if crit==0
  crit=1e-6;
end

I = DimX(1);
J = DimX(2);
K = DimX(3);

InitWithRandom=0;
if nargin<8
   InitWithRandom=1;
end
if nargin>7 & size(A,1)~=I
  InitWithRandom=1;
end

if nargin<5
   ConstA = 0;ConstB = 0;ConstC = 0;
else
   ConstA = Constraints(1);ConstB = Constraints(2);ConstC = Constraints(3);
end

if InitWithRandom

  if I<Fac
    A = rand(I,Fac);
  else
    A = orth(rand(I,Fac));
  end
  if J<Fac
    B = rand(J,Fac);
  else
    B = orth(rand(J,Fac));
  end
  if K<Fac
    C = rand(K,Fac);
  else
    C = orth(rand(K,Fac));
  end
end

SumSqX = sum(sum(abs(X).^2));
fit    = SumSqX;
fit0   = fit;
fitold = 2*fit;
it     = 0;
Delta  = 5;

while abs((fit-fitold)/fitold)>crit&it<maxit&fit>10*eps
   it=it+1;
   fitold=fit;

   % Do line-search
   if rem(it+2,2)==-1
      [A,B,C,Delta]=linesrch(X,DimX,A,B,C,Ao,Bo,Co,Delta);
   end
   
   Ao=A;Bo=B;Co=C;
   % Update A
   Xbc=0;
   for k=1:K
     Xbc = Xbc + X(:,(k-1)*J+1:k*J)*conj(B*diag(C(k,:)));
   end
   if ConstA == 0 % Unconstrained
      A = Xbc*pinv((B'*B).*(C'*C)).';
   elseif ConstA == 1 % Nonnegativity, requires reals
      Aold = A;
      for i = 1:I
         ztz = (B'*B).*(C'*C);
         A(i,:) = fastnnls(ztz,Xbc(i,:)')';
      end
      if any(sum(A)<100*eps*I)
         A = .99*Aold+.01*A; % To prevent a matrix with zero columns
      end
   elseif ConstA == 2 % Orthogonality
      A = Xbc*(Xbc'*Xbc)^(-.5);
   elseif ConstA == 3 % Unimodality
      A = unimodalcrossproducts((B'*B).*(C'*C),Xbc',A);
   end

   % Project X down on orth(A) - saves time if first mode is large
   [Qa,Ra]=qr(A,0);
   x=Qa'*X;

   % Update B
   if ConstB == 10 % Procrustes
      B = eye(Fac);
   else
      Xac=0;
      for k=1:K
         Xac = Xac + x(:,(k-1)*J+1:k*J).'*conj(Ra*diag(C(k,:)));
      end
      if ConstB == 0 % Unconstrained
         B = Xac*pinv((Ra'*Ra).*(C'*C)).';
      elseif ConstB == 1 % Nonnegativity, requires reals
         Bold = B;
         for j = 1:J
            ztz = (Ra'*Ra).*(C'*C);
            B(j,:) = fastnnls(ztz,Xac(j,:)')';
         end
         if any(sum(B)<100*eps*J)
            B = .99*Bold+.01*B; % To prevent a matrix with zero columns
         end
      end
   end
  
    % Update C
    if ConstC == 0 % Unconstrained
       ab=pinv((Ra'*Ra).*(B'*B));
       for k=1:K 
          C(k,:) = (ab*diag(Ra'* x(:,(k-1)*J+1:k*J)*conj(B))).';
       end
    elseif ConstC == 1  % Nonnegativity, requires reals
       Cold = C;
       ztz = (Ra'*Ra).*(B'*B);
       for k = 1:K
          xab = diag(Ra'* x(:,(k-1)*J+1:k*J)*B);
          C(k,:) = fastnnls(ztz,xab)';
       end
       if any(sum(C)<100*eps*K)
          C = .99*Cold+.01*C; % To prevent a matrix with zero columns
       end
    elseif ConstC == 2 % Orthogonality
       Z=(Ra'*Ra).*(B'*B);
       Y=[];
       for k=1:K
          d=diag(Ra'*x(:,(k-1)*J+1:k*J)*B)'; 
          Y=[Y;d];
       end;
       [P,D,Q]=svd(Y,0);
       C=P*Q';
    elseif ConstC == 3 % Unimodality
       xab = [];
       for k = 1:K
          xab = [xab diag(Ra'* x(:,(k-1)*J+1:k*J)*B)];
       end
       C = unimodalcrossproducts((Ra'*Ra).*(B'*B),xab,C);
    elseif ConstC == 10 % GPA => Isotropic scaling factor
       ab=(Ra'*Ra).*(B'*B);
       ab = pinv(ab(:));
       C(1,:) = 1;
       for k=2:K 
          yy = [];
          yyy = diag(Ra'* x(:,(k-1)*J+1:k*J)*conj(B)).';
          for f=1:Fac
             yy = [yy;yyy(:)];
          end
          C(k,:) = ab*yy;
       end
    end
      
    % Calculating fit. Using orthogonalization instead
   %fit=0;for k=1:K,residual=X(:,(k-1)*J+1:k*J)-A*diag(C(k,:))*B.';fit=fit+sum(sum((abs(residual).^2)));end
   [Qb,Rb]=qr(B,0);
   [Z,Rc]=qr(C,0);
   fit=SumSqX-sum(sum(abs(Ra*ppp(Rb,Rc).').^2));
   
   if rem(it,showfit)==0
      fprintf(' %12.10f       %g        %3.4f \n',fit,it,100*(1-fit/fit0));
   end
end

% ORDER ACCORDING TO VARIANCE
Tuck     = diag((A'*A).*(B'*B).*(C'*C));
[out,ID] = sort(Tuck);
A        = A(:,ID);
if ConstB ~= 10 % Else B is eye
   B        = B(:,ID);
end
C        = C(:,ID);
% NORMALIZE A AND C (variance in B)
if ConstB ~= 10 % Then B is eye
   for f=1:Fac,normC(f) = norm(C(:,f));end
   for f=1:Fac,normA(f) = norm(A(:,f));end
   B        = B*diag(normC)*diag(normA);  
   A        = A*diag(normA.^(-1));
   C        = C*diag(normC.^(-1));
   
   % APPLY SIGN CONVENTION
   SignA = sign(sum(sign(A))+eps);
   SignC = sign(sum(sign(C))+eps);
   A = A*diag(SignA);
   C = C*diag(SignC);
   B = B*diag(SignA)*diag(SignC);
end

function [NewA,NewB,NewC,DeltaMin] = linesrch(X,DimX,A,B,C,Ao,Bo,Co,Delta);

dbg=0;

if nargin<5
  Delta=5;
else
  Delta=max(2,Delta);
end

dA=A-Ao;
dB=B-Bo;
dC=C-Co;
Fit1=sum(sum(abs(X-A*ppp(B,C).').^2));
regx=[1 0 0 Fit1];
Fit2=sum(sum(abs(X-(A+Delta*dA)*ppp((B+Delta*dB),(C+Delta*dC)).').^2));
regx=[regx;1 Delta Delta.^2 Fit2];

while Fit2>Fit1
  if dbg
    disp('while Fit2>Fit1')
  end
  Delta=Delta*.6;
  Fit2=sum(sum(abs(X-(A+Delta*dA)*ppp((B+Delta*dB),(C+Delta*dC)).').^2));
  regx=[regx;1 Delta Delta.^2 Fit2];
end

Fit3=sum(sum(abs(X-(A+2*Delta*dA)*ppp((B+2*Delta*dB),(C+2*Delta*dC)).').^2));
regx=[regx;1 2*Delta (2*Delta).^2 Fit3];

while Fit3<Fit2
  if dbg
    disp('while Fit3<Fit2')
  end
  Delta=1.8*Delta;
  Fit2=Fit3;
  Fit3=sum(sum(abs(X-(A+2*Delta*dA)*ppp((B+2*Delta*dB),(C+2*Delta*dC)).').^2));
  regx=[regx;1 2*Delta (2*Delta).^2 Fit2];
end

% Add one point between the two smallest fits
[a,b]=sort(regx(:,4));
regx=regx(b,:);
Delta4=(regx(1,2)+regx(2,2))/2;
Fit4=sum(sum(abs(X-(A+Delta4*dA)*ppp((B+Delta4*dB),(C+Delta4*dC)).').^2));
regx=[regx;1 Delta4 Delta4.^2 Fit4];

%reg=pinv([1 0 0;1 Delta Delta^2;1 2*Delta (2*Delta)^2])*[Fit1;Fit2;Fit3]
reg=pinv(regx(:,1:3))*regx(:,4);
%DeltaMin=2*reg(3);

DeltaMin=-reg(2)/(2*reg(3));

%a*x2 + bx + c = fit
%2ax + b = 0
%x=-b/2a

NewA=A+DeltaMin*dA;
NewB=B+DeltaMin*dB;
NewC=C+DeltaMin*dC;
Fit=sum(sum(abs(X-NewA*ppp(NewB,NewC).').^2));

if dbg
  regx
  plot(regx(:,2),regx(:,4),'o'),
  hold on
  x=linspace(0,max(regx(:,2))*1.2);
  plot(x',[ones(100,1) x' x'.^2]*reg),
  hold off
  drawnow
  [DeltaMin Fit],pause
end

[minfit,number]=min(regx(:,4));
if Fit>minfit
  DeltaMin=regx(number,2);
  NewA=A+DeltaMin*dA;
  NewB=B+DeltaMin*dB;
  NewC=C+DeltaMin*dC;
end

function AB=ppp(A,B);

% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $
%
% Copyright, 1998 - 
% This M-file and the code in it belongs to the holder of the
% copyrights and is made public under the following constraints:
% It must not be changed or modified and code cannot be added.
% The file must be regarded as read-only. Furthermore, the
% code can not be made part of anything but the 'N-way Toolbox'.
% In case of doubt, contact the holder of the copyrights.
%
% Rasmus Bro
% Chemometrics Group, Food Technology
% Department of Food and Dairy Science
% Royal Veterinary and Agricultutal University
% Rolighedsvej 30, DK-1958 Frederiksberg, Denmark
% Phone  +45 35283296
% Fax    +45 35283245
% E-mail rb@kvl.dk
%
% The parallel proportional profiles product - triple-P product
% For two matrices with similar column dimension the triple-P product
% is ppp(A,B) = [kron(B(:,1),A(:,1) .... kron(B(:,F),A(:,F)]
% 
% AB = ppp(A,B);
%
% Copyright 1998
% Rasmus Bro
% KVL,DK
% rb@kvl.dk

[I,F]=size(A);
[J,F1]=size(B);

if F~=F1
   error(' Error in ppp.m - The matrices must have the same number of columns')
end

AB=zeros(I*J,F);
for f=1:F
   ab=A(:,f)*B(:,f).';
   AB(:,f)=ab(:);
end



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
%	See also LSCOV, SLASH.

%	L. Shure 5-8-87
%	Revised, 12-15-88,8-31-89 LS.
%	Copyright (c) 1984-94 by The MathWorks, Inc.

%       Revised by:
%	Copyright
%	Rasmus Bro 1995
%	Denmark
%	E-mail rb@kvl.dk
%       According to Bro & de Jong, J. Chemom, 1997

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


function B=unimodalcrossproducts(XtX,XtY,Bold)

% Solves the problem min|Y-XB'| subject to the columns of 
% B are unimodal and nonnegative. The algorithm is iterative and
% only one iteration is given, hence the solution is only improving 
% the current estimate
%
% I/O B=unimodalcrossproducts(XtX,XtY,Bold)
% Modified from unimodal.m to handle crossproducts in input 1999
%
% Copyright 1997
%
% Rasmus Bro
% Royal Veterinary & Agricultural University
% Denmark
% rb@kvl.dk
%
% Reference
% Bro and Sidiropoulos, "Journal of Chemometrics", 1998, 12, 223-247. 


B=Bold;
F=size(B,2);
for f=1:F
   xty = XtY(f,:)-XtX(f,[1:f-1 f+1:F])*B(:,[1:f-1 f+1:F])';
   beta=pinv(XtX(f,f))*xty;
   B(:,f)=ulsr(beta',1);
end


function [b,All,MaxML]=ulsr(x,NonNeg);

% ------INPUT------
%
% x          is the vector to be approximated
% NonNeg     If NonNeg is one, nonnegativity is imposed
%
%
%
% ------OUTPUT-----
%
% b 	     is the best ULSR vector
% All 	     is containing in its i'th column the ULSRFIX solution for mode
% 	     location at the i'th element. The ULSR solution given in All
%            is found disregarding the i'th element and hence NOT optimal
% MaxML      is the optimal (leftmost) mode location (i.e. position of maximum)
%
% ___________________________________________________________
%
%
%               Copyright 1997
%
% Nikos Sidiroupolos
% University of Maryland
% Maryland, US
%
%       &
%
% Rasmus Bro
% Royal Veterinary & Agricultural University
% Denmark
%
% 
% ___________________________________________________________


% This file uses MONREG.M

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

function [b,B,AllBs]=monreg(x);

% Monotonic regression according
% to J. B. Kruskal 64
%
% b     = min|x-b| subject to monotonic increase
% B     = b, but condensed
% AllBs = All monotonic regressions, i.e. AllBs(1:i,i) is the 
%         monotonic regression of x(1:i)
%
%
% Copyright 1997
%
% Rasmus Bro
% Royal Veterinary & Agricultural University
% Denmark
% rb@kvl.dk
%


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
