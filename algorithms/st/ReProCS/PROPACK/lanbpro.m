function [U,B_k,V,p,ierr,work] = lanbpro(varargin)

%LANBPRO Lanczos bidiagonalization with partial reorthogonalization.
%   LANBPRO computes the Lanczos bidiagonalization of a real 
%   matrix using the  with partial reorthogonalization. 
%
%   [U_k,B_k,V_k,R,ierr,work] = LANBPRO(A,K,R0,OPTIONS,U_old,B_old,V_old) 
%   [U_k,B_k,V_k,R,ierr,work] = LANBPRO('Afun','Atransfun',M,N,K,R0, ...
%                                       OPTIONS,U_old,B_old,V_old) 
%
%   Computes K steps of the Lanczos bidiagonalization algorithm with partial 
%   reorthogonalization (BPRO) with M-by-1 starting vector R0, producing a 
%   lower bidiagonal K-by-K matrix B_k, an N-by-K matrix V_k, an M-by-K 
%   matrix U_k and an M-by-1 vector R such that
%        A*V_k = U_k*B_k + R
%   Partial reorthogonalization is used to keep the columns of V_K and U_k
%   semiorthogonal:
%         MAX(DIAG((EYE(K) - V_K'*V_K))) <= OPTIONS.delta 
%   and 
%         MAX(DIAG((EYE(K) - U_K'*U_K))) <= OPTIONS.delta.
%
%   B_k = LANBPRO(...) returns the bidiagonal matrix only.
%
%   The first input argument is either a real matrix, or a string
%   containing the name of an M-file which applies a linear operator 
%   to the columns of a given matrix. In the latter case, the second 
%   input must be the name of an M-file which applies the transpose of 
%   the same linear operator to the columns of a given matrix,  
%   and the third and fourth arguments must be M and N, the dimensions 
%   of then problem.
%
%   The OPTIONS structure is used to control the reorthogonalization:
%     OPTIONS.delta:  Desired level of orthogonality 
%                     (default = sqrt(eps/K)).
%     OPTIONS.eta  :  Level of orthogonality after reorthogonalization 
%                     (default = eps^(3/4)/sqrt(K)).
%     OPTIONS.cgs  :  Flag for switching between different reorthogonalization
%                     algorithms:
%                      0 = iterated modified Gram-Schmidt  (default)
%                      1 = iterated classical Gram-Schmidt 
%     OPTIONS.elr  :  If OPTIONS.elr = 1 (default) then extended local
%                     reorthogonalization is enforced.
%     OPTIONS.onesided
%                  :  If OPTIONS.onesided = 0 (default) then both the left
%                     (U) and right (V) Lanczos vectors are kept 
%                     semiorthogonal. 
%                     OPTIONS.onesided = 1 then only the columns of U are
%                     are reorthogonalized.
%                     OPTIONS.onesided = -1 then only the columns of V are
%                     are reorthogonalized.
%     OPTIONS.waitbar
%                  :  The progress of the algorithm is display graphically.
%
%   If both R0, U_old, B_old, and V_old are provided, they must
%   contain a partial Lanczos bidiagonalization of A on the form
%
%        A V_old = U_old B_old + R0 .  
%
%   In this case the factorization is extended to dimension K x K by
%   continuing the Lanczos bidiagonalization algorithm with R0 as a 
%   starting vector.
%
%   The output array work contains information about the work used in
%   reorthogonalizing the u- and v-vectors.
%      work = [ RU  PU ]
%             [ RV  PV ] 
%   where
%      RU = Number of reorthogonalizations of U.
%      PU = Number of inner products used in reorthogonalizing U.
%      RV = Number of reorthogonalizations of V.
%      PV = Number of inner products used in reorthogonalizing V.

% References: 
% R.M. Larsen, Ph.D. Thesis, Aarhus University, 1998.
%
% G. H. Golub & C. F. Van Loan, "Matrix Computations",
% 3. Ed., Johns Hopkins, 1996.  Section 9.3.4.
%
% B. N. Parlett, ``The Symmetric Eigenvalue Problem'', 
% Prentice-Hall, Englewood Cliffs, NJ, 1980.
%
% H. D. Simon, ``The Lanczos algorithm with partial reorthogonalization'',
% Math. Comp. 42 (1984), no. 165, 115--142.
%

% Rasmus Munk Larsen, DAIMI, 1998.

% Check input arguments.

global LANBPRO_TRUTH
LANBPRO_TRUTH=0;

if LANBPRO_TRUTH==1
  global MU NU MUTRUE NUTRUE
  global MU_AFTER NU_AFTER MUTRUE_AFTER NUTRUE_AFTER
end

if nargin<1 | length(varargin)<2
  error('Not enough input arguments.');
end
narg=length(varargin);

A = varargin{1};
if isnumeric(A) | isstruct(A)
  if isnumeric(A)
    if ~isreal(A)
      error('A must be real')
    end  
    [m n] = size(A);
  elseif isstruct(A)
    [m n] = size(A.R);
  end
  k=varargin{2};
  if narg >= 3 & ~isempty(varargin{3});
    p = varargin{3};
  else
    p = rand(m,1)-0.5;
  end
  if narg < 4, options = []; else options=varargin{4}; end
  if narg > 4 
    if narg<7
      error('All or none of U_old, B_old and V_old must be provided.')
    else
      U = varargin{5}; B_k = varargin{6}; V = varargin{7};
    end
  else
    U = []; B_k = []; V = [];
  end
  if narg > 7, anorm=varargin{8}; else anorm = []; end
else
  if narg<5
    error('Not enough input arguments.');
  end
  Atrans = varargin{2};
  if ~isstr(Atrans)
    error('Afunc and Atransfunc must be names of m-files')
  end
  m = varargin{3};
  n = varargin{4};
  if ~isreal(n) | abs(fix(n)) ~= n | ~isreal(m) | abs(fix(m)) ~= m
    error('M and N must be positive integers.')
  end
  k=varargin{5};
  if narg < 6, p = rand(m,1)-0.5; else p=varargin{6}; end  
  if narg < 7, options = []; else options=varargin{7}; end  
  if narg > 7
    if  narg < 10
      error('All or none of U_old, B_old and V_old must be provided.')
    else
      U = varargin{8}; B_k = varargin{9}; V = varargin{10};
    end
  else
    U = []; B_k = []; V=[];
  end
  if narg > 10, anorm=varargin{11}; else anorm = [];  end
end

% Quick return for min(m,n) equal to 0 or 1.
if min(m,n) == 0
   U = [];  B_k = [];  V = [];  p = [];  ierr = 0;  work = zeros(2,2);
   return
elseif  min(m,n) == 1
  if isnumeric(A)
    U = 1;  B_k = A;  V = 1;  p = 0; ierr = 0; work = zeros(2,2);
  else
    U = 1;  B_k = feval(A,1); V = 1; p = 0; ierr = 0; work = zeros(2,2);
  end
  if nargout<3
    U = B_k;
  end
  return
end

% Set options.  
%m2 = 3/2*(sqrt(m)+1);
%n2 = 3/2*(sqrt(n)+1);
m2 = 3/2;
n2 = 3/2;
delta = sqrt(eps/k); % Desired level of orthogonality.
eta = eps^(3/4)/sqrt(k);    % Level of orth. after reorthogonalization.
cgs = 0;             % Flag for switching between iterated MGS and CGS.
elr = 2;             % Flag for switching extended local 
                     % reorthogonalization on and off.
gamma = 1/sqrt(2);   % Tolerance for iterated Gram-Schmidt.
onesided = 0; t = 0; waitb = 0;

% Parse options struct
if ~isempty(options) & isstruct(options)
  c = fieldnames(options);
  for i=1:length(c)
    if strmatch(c(i),'delta'), delta = getfield(options,'delta');  end
    if strmatch(c(i),'eta'), eta = getfield(options,'eta'); end
    if strmatch(c(i),'cgs'), cgs = getfield(options,'cgs'); end
    if strmatch(c(i),'elr'), elr = getfield(options,'elr'); end
    if strmatch(c(i),'gamma'), gamma = getfield(options,'gamma'); end
    if strmatch(c(i),'onesided'), onesided = getfield(options,'onesided'); end
    if strmatch(c(i),'waitbar'), waitb=1; end
  end
end

if waitb
  waitbarh = waitbar(0,'Lanczos bidiagonalization in progress...');
end

if isempty(anorm)
  anorm = []; est_anorm=1; 
else
  est_anorm=0; 
end

% Conservative statistical estimate on the size of round-off terms. 
% Notice that {\bf u} == eps/2.
FUDGE = 1.01; % Fudge factor for ||A||_2 estimate.

npu = 0; npv = 0; ierr = 0;
p = p(:);
% Prepare for Lanczos iteration.
if isempty(U)
  V = zeros(n,k); U = zeros(m,k);
  beta = zeros(k+1,1); alpha = zeros(k,1);
  beta(1) = norm(p);
  % Initialize MU/NU-recurrences for monitoring loss of orthogonality.
  nu = zeros(k,1); mu = zeros(k+1,1);
  mu(1)=1; nu(1)=1;
  
  numax = zeros(k,1); mumax = zeros(k,1);
  force_reorth = 0;  nreorthu = 0; nreorthv = 0;
  j0 = 1;
else
  j = size(U,2); % Size of existing factorization
  % Allocate space for Lanczos vectors
  U = [U, zeros(m,k-j)];
  V = [V, zeros(n,k-j)];
  alpha = zeros(k+1,1);  beta = zeros(k+1,1);
  alpha(1:j) = diag(B_k); if j>1 beta(2:j) = diag(B_k,-1); end
  beta(j+1) = norm(p);
  % Reorthogonalize p.
  if j<k & beta(j+1)*delta < anorm*eps,
    fro = 1;
    ierr = j;
  end
  int = [1:j]';
  [p,beta(j+1),rr] = reorth(U,p,beta(j+1),int,gamma,cgs);
  npu =  rr*j;  nreorthu = 1;  force_reorth= 1;  

  % Compute Gerscgorin bound on ||B_k||_2
  if est_anorm
    anorm = FUDGE*sqrt(norm(B_k'*B_k,1));
  end
  mu = m2*eps*ones(k+1,1); nu = zeros(k,1);
  numax = zeros(k,1); mumax = zeros(k,1);
  force_reorth = 1;  nreorthu = 0; nreorthv = 0;
  j0 = j+1;
end


if isnumeric(A)
  At = A';
end

if delta==0
  fro = 1; % The user has requested full reorthogonalization.
else
  fro = 0;
end

if LANBPRO_TRUTH==1
  MUTRUE = zeros(k,k); NUTRUE = zeros(k-1,k-1);
  MU = zeros(k,k); NU = zeros(k-1,k-1);
  
  MUTRUE_AFTER = zeros(k,k); NUTRUE_AFTER = zeros(k-1,k-1);
  MU_AFTER = zeros(k,k); NU_AFTER = zeros(k-1,k-1);
end

% Perform Lanczos bidiagonalization with partial reorthogonalization.
for j=j0:k
  if waitb
    waitbar(j/k,waitbarh)
  end

  if beta(j) ~= 0
    U(:,j) = p/beta(j);
  else
    U(:,j) = p;
  end

  % Replace norm estimate with largest Ritz value.
  if j==6
    B = [[diag(alpha(1:j-1))+diag(beta(2:j-1),-1)]; ...
      [zeros(1,j-2),beta(j)]];
    anorm = FUDGE*norm(B);
    est_anorm = 0;
  end
  
  %%%%%%%%%% Lanczos step to generate v_j. %%%%%%%%%%%%%
  if j==1
    if isnumeric(A)
      r = At*U(:,1);    
    elseif isstruct(A)
      r = A.R\U(:,1);          
    else
      r = feval(Atrans,U(:,1));
    end
    alpha(1) = norm(r);
    if est_anorm
      anorm = FUDGE*alpha(1);
    end
  else    
    if isnumeric(A)
      r = At*U(:,j) - beta(j)*V(:,j-1);
    elseif isstruct(A)
      r = A.R\U(:,j) - beta(j)*V(:,j-1);      
    else
      r = feval(Atrans,U(:,j))  - beta(j)*V(:,j-1);
    end
    alpha(j) = norm(r); 

    % Extended local reorthogonalization    
    if alpha(j)<gamma*beta(j) & elr & ~fro
      normold = alpha(j);
      stop = 0;
      while ~stop
	t = V(:,j-1)'*r;
	r = r - V(:,j-1)*t;
	alpha(j) = norm(r);
	if beta(j) ~= 0
	  beta(j) = beta(j) + t;
	end
	if alpha(j)>=gamma*normold
	  stop = 1;
	else
	  normold = alpha(j);
	end
      end
    end

    if est_anorm
      if j==2
	anorm = max(anorm,FUDGE*sqrt(alpha(1)^2+beta(2)^2+alpha(2)*beta(2)));
      else	
	anorm = max(anorm,FUDGE*sqrt(alpha(j-1)^2+beta(j)^2+alpha(j-1)* ...
	    beta(j-1) + alpha(j)*beta(j)));
      end			     
    end
    
    if ~fro & alpha(j) ~= 0
      % Update estimates of the level of orthogonality for the
      %  columns 1 through j-1 in V.
      nu = update_nu(nu,mu,j,alpha,beta,anorm);
      numax(j) = max(abs(nu(1:j-1)));
    end

    if j>1 & LANBPRO_TRUTH
      NU(1:j-1,j-1) = nu(1:j-1);
      NUTRUE(1:j-1,j-1) = V(:,1:j-1)'*r/alpha(j);
    end
    
    if elr>0
      nu(j-1) = n2*eps;
    end
    
    % IF level of orthogonality is worse than delta THEN 
    %    Reorthogonalize v_j against some previous  v_i's, 0<=i<j.
    if onesided~=-1 & ( fro | numax(j) > delta | force_reorth ) & alpha(j)~=0
      % Decide which vectors to orthogonalize against:
      if fro | eta==0
	int = [1:j-1]';
      elseif force_reorth==0
	int = compute_int(nu,j-1,delta,eta,0,0,0);
      end
      % Else use int from last reorth. to avoid spillover from mu_{j-1} 
      % to nu_j.
      
      % Reorthogonalize v_j 
      [r,alpha(j),rr] = reorth(V,r,alpha(j),int,gamma,cgs);
      npv = npv + rr*length(int); % number of inner products.
      nu(int) = n2*eps;  % Reset nu for orthogonalized vectors.

      % If necessary force reorthogonalization of u_{j+1} 
      % to avoid spillover
      if force_reorth==0 
	force_reorth = 1; 
      else
	force_reorth = 0; 
      end
      nreorthv = nreorthv + 1;
    end
  end

  
  % Check for convergence or failure to maintain semiorthogonality
  if alpha(j) < max(n,m)*anorm*eps & j<k, 
    % If alpha is "small" we deflate by setting it
    % to 0 and attempt to restart with a basis for a new 
    % invariant subspace by replacing r with a random starting vector:
    %j
    %disp('restarting, alpha = 0')
    alpha(j) = 0;
    bailout = 1;
    for attempt=1:3    
      r = rand(m,1)-0.5;  
      if isnumeric(A)
	r = At*r;    
      elseif isstruct(A)
	r = A.R\r;    
      else
	r = feval(Atrans,r);
      end
      nrm=sqrt(r'*r); % not necessary to compute the norm accurately here.
      int = [1:j-1]';
      [r,nrmnew,rr] = reorth(V,r,nrm,int,gamma,cgs);
      npv = npv + rr*length(int(:));        nreorthv = nreorthv + 1;
      nu(int) = n2*eps;
      if nrmnew > 0
	% A vector numerically orthogonal to span(Q_k(:,1:j)) was found. 
	% Continue iteration.
	bailout=0;
	break;
      end
    end
    if bailout
      j = j-1;
      ierr = -j;
      break;
    else
      r=r/nrmnew; % Continue with new normalized r as starting vector.
      force_reorth = 1;
      if delta>0
	fro = 0;    % Turn off full reorthogonalization.
      end
    end       
  elseif  j<k & ~fro & anorm*eps > delta*alpha(j)
%    fro = 1;
    ierr = j;
  end

  if j>1 & LANBPRO_TRUTH
    NU_AFTER(1:j-1,j-1) = nu(1:j-1);
    NUTRUE_AFTER(1:j-1,j-1) = V(:,1:j-1)'*r/alpha(j);
  end

  
  if alpha(j) ~= 0
    V(:,j) = r/alpha(j);
  else
    V(:,j) = r;
  end

  %%%%%%%%%% Lanczos step to generate u_{j+1}. %%%%%%%%%%%%%
  if waitb
    waitbar((2*j+1)/(2*k),waitbarh)
  end
  
  if isnumeric(A)
    p = A*V(:,j) - alpha(j)*U(:,j);
  elseif isstruct(A)
    p = A.Rt\V(:,j) - alpha(j)*U(:,j);
  else
    p = feval(A,V(:,j)) - alpha(j)*U(:,j);
  end
  beta(j+1) = norm(p);
  % Extended local reorthogonalization
  if beta(j+1)<gamma*alpha(j) & elr & ~fro
    normold = beta(j+1);
    stop = 0;
    while ~stop
      t = U(:,j)'*p;
      p = p - U(:,j)*t;
      beta(j+1) = norm(p);
      if alpha(j) ~= 0 
	alpha(j) = alpha(j) + t;
      end
      if beta(j+1) >= gamma*normold
	stop = 1;
      else
	normold = beta(j+1);
      end
    end
  end

  if est_anorm
    % We should update estimate of ||A||  before updating mu - especially  
    % important in the first step for problems with large norm since alpha(1) 
    % may be a severe underestimate!  
    if j==1
      anorm = max(anorm,FUDGE*pythag(alpha(1),beta(2))); 
    else
      anorm = max(anorm,FUDGE*sqrt(alpha(j)^2+beta(j+1)^2 + alpha(j)*beta(j)));
    end
  end
  
  
  if ~fro & beta(j+1) ~= 0
    % Update estimates of the level of orthogonality for the columns of V.
    mu = update_mu(mu,nu,j,alpha,beta,anorm);
    mumax(j) = max(abs(mu(1:j)));  
  end

  if LANBPRO_TRUTH==1
    MU(1:j,j) = mu(1:j);
    MUTRUE(1:j,j) = U(:,1:j)'*p/beta(j+1);
  end
  
  if elr>0
    mu(j) = m2*eps;
  end
  
  % IF level of orthogonality is worse than delta THEN 
  %    Reorthogonalize u_{j+1} against some previous  u_i's, 0<=i<=j.
  if onesided~=1 & (fro | mumax(j) > delta | force_reorth) & beta(j+1)~=0
    % Decide which vectors to orthogonalize against.
    if fro | eta==0
      int = [1:j]';
    elseif force_reorth==0
      int = compute_int(mu,j,delta,eta,0,0,0); 
    else
      int = [int; max(int)+1];
    end
    % Else use int from last reorth. to avoid spillover from nu to mu.

%    if onesided~=0
%      fprintf('i = %i, nr = %i, fro = %i\n',j,size(int(:),1),fro)
%    end
    % Reorthogonalize u_{j+1} 
    [p,beta(j+1),rr] = reorth(U,p,beta(j+1),int,gamma,cgs);    
    npu = npu + rr*length(int);  nreorthu = nreorthu + 1;    

    % Reset mu to epsilon.
    mu(int) = m2*eps;    
    
    if force_reorth==0 
      force_reorth = 1; % Force reorthogonalization of v_{j+1}.
    else
      force_reorth = 0; 
    end
  end
  
  % Check for convergence or failure to maintain semiorthogonality
  if beta(j+1) < max(m,n)*anorm*eps  & j<k,     
    % If beta is "small" we deflate by setting it
    % to 0 and attempt to restart with a basis for a new 
    % invariant subspace by replacing p with a random starting vector:
    %j
    %disp('restarting, beta = 0')
    beta(j+1) = 0;
    bailout = 1;
    for attempt=1:3    
      p = rand(n,1)-0.5;  
      if isnumeric(A)
	p = A*p;
      elseif isstruct(A)
	p = A.Rt\p;
      else
	p = feval(A,p);
      end
      nrm=sqrt(p'*p); % not necessary to compute the norm accurately here.
      int = [1:j]';
      [p,nrmnew,rr] = reorth(U,p,nrm,int,gamma,cgs);
      npu = npu + rr*length(int(:));  nreorthu = nreorthu + 1;
      mu(int) = m2*eps;
      if nrmnew > 0
	% A vector numerically orthogonal to span(Q_k(:,1:j)) was found. 
	% Continue iteration.
	bailout=0;
	break;
      end
    end
    if bailout
      ierr = -j;
      break;
    else
      p=p/nrmnew; % Continue with new normalized p as starting vector.
      force_reorth = 1;
      if delta>0
	fro = 0;    % Turn off full reorthogonalization.
      end
    end       
  elseif  j<k & ~fro & anorm*eps > delta*beta(j+1) 
%    fro = 1;
    ierr = j;
  end  
  
  if LANBPRO_TRUTH==1
    MU_AFTER(1:j,j) = mu(1:j);
    MUTRUE_AFTER(1:j,j) = U(:,1:j)'*p/beta(j+1);
  end  
end
if waitb
  close(waitbarh)
end

if j<k
  k = j;
end

B_k = spdiags([alpha(1:k) [beta(2:k);0]],[0 -1],k,k);
if nargout==1
  U = B_k;
elseif k~=size(U,2) | k~=size(V,2)  
  U = U(:,1:k);
  V = V(:,1:k);
end
if nargout>5
  work = [[nreorthu,npu];[nreorthv,npv]];
end



function mu = update_mu(muold,nu,j,alpha,beta,anorm)

% UPDATE_MU:  Update the mu-recurrence for the u-vectors.
%
%   mu_new = update_mu(mu,nu,j,alpha,beta,anorm)

%  Rasmus Munk Larsen, DAIMI, 1998.

binv = 1/beta(j+1);
mu = muold;
eps1 = 100*eps/2;
if j==1
  T = eps1*(pythag(alpha(1),beta(2)) + pythag(alpha(1),beta(1)));
  T = T + eps1*anorm;
  mu(1) = T / beta(2);
else
  mu(1) = alpha(1)*nu(1) - alpha(j)*mu(1);
%  T = eps1*(pythag(alpha(j),beta(j+1)) + pythag(alpha(1),beta(1)));
  T = eps1*(sqrt(alpha(j).^2+beta(j+1).^2) + sqrt(alpha(1).^2+beta(1).^2));
  T = T + eps1*anorm;
  mu(1) = (mu(1) + sign(mu(1))*T) / beta(j+1);
  % Vectorized version of loop:
  if j>2
    k=2:j-1;
    mu(k) = alpha(k).*nu(k) + beta(k).*nu(k-1) - alpha(j)*mu(k);
    %T = eps1*(pythag(alpha(j),beta(j+1)) + pythag(alpha(k),beta(k)));
    T = eps1*(sqrt(alpha(j).^2+beta(j+1).^2) + sqrt(alpha(k).^2+beta(k).^2));
    T = T + eps1*anorm;
    mu(k) = binv*(mu(k) + sign(mu(k)).*T);
  end
%  T = eps1*(pythag(alpha(j),beta(j+1)) + pythag(alpha(j),beta(j)));
  T = eps1*(sqrt(alpha(j).^2+beta(j+1).^2) + sqrt(alpha(j).^2+beta(j).^2));
  T = T + eps1*anorm;
  mu(j) = beta(j)*nu(j-1);
  mu(j) = (mu(j) + sign(mu(j))*T) / beta(j+1);
end  
mu(j+1) = 1;


function nu = update_nu(nuold,mu,j,alpha,beta,anorm)

% UPDATE_MU:  Update the nu-recurrence for the v-vectors.
%
%  nu_new = update_nu(nu,mu,j,alpha,beta,anorm)

%  Rasmus Munk Larsen, DAIMI, 1998.

nu = nuold;
ainv = 1/alpha(j);
eps1 = 100*eps/2;
if j>1
  k = 1:(j-1);
%  T = eps1*(pythag(alpha(k),beta(k+1)) + pythag(alpha(j),beta(j)));
  T = eps1*(sqrt(alpha(k).^2+beta(k+1).^2) + sqrt(alpha(j).^2+beta(j).^2));
  T = T + eps1*anorm;
  nu(k) = beta(k+1).*mu(k+1) + alpha(k).*mu(k) - beta(j)*nu(k);
  nu(k) = ainv*(nu(k) + sign(nu(k)).*T);
end
nu(j) = 1;

function x = pythag(y,z)
%PYTHAG Computes sqrt( y^2 + z^2 ).
%
% x = pythag(y,z)
%
% Returns sqrt(y^2 + z^2) but is careful to scale to avoid overflow.

% Christian H. Bischof, Argonne National Laboratory, 03/31/89.

[m n] = size(y);
if m>1 | n>1
  y = y(:); z=z(:);
  rmax = max(abs([y z]'))';
  id=find(rmax==0);
  if length(id)>0
    rmax(id) = 1;
    x = rmax.*sqrt((y./rmax).^2 + (z./rmax).^2);
    x(id)=0;
  else
    x = rmax.*sqrt((y./rmax).^2 + (z./rmax).^2);
  end
  x = reshape(x,m,n);
else
  rmax = max(abs([y;z]));
  if (rmax==0)
    x = 0;
  else
    x = rmax*sqrt((y/rmax)^2 + (z/rmax)^2);
  end
end
  
