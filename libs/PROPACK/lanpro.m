function [Q_k,T_k,r,anorm,ierr,work] = lanpro(A,nin,kmax,r,options,...
    Q_k,T_k,anorm)
 
%LANPRO   Lanczos tridiagonalization with partial reorthogonalization
%   LANPRO computes the Lanczos tridiagonalization of a real symmetric 
%   matrix using the symmetric Lanczos algorithm with partial 
%   reorthogonalization. 
%
%   [Q_K,T_K,R,ANORM,IERR,WORK] = LANPRO(A,K,R0,OPTIONS,Q_old,T_old)
%   [Q_K,T_K,R,ANORM,IERR,WORK] = LANPRO('Afun',N,K,R0,OPTIONS,Q_old,T_old)
%
%   Computes K steps of the Lanczos algorithm with starting vector R0, 
%   and returns the K x K tridiagonal T_K, the N x K matrix Q_K 
%   with semiorthonormal columns and the residual vector R such that 
%
%        A*Q_K = Q_K*T_K + R .
%
%   Partial reorthogonalization is used to keep the columns of Q_K 
%   semiorthogonal:
%        MAX(DIAG((eye(k) - Q_K'*Q_K))) <= OPTIONS.delta.
%
%
%   The first input argument is either a real symmetric matrix, a struct with
%   components A.L and A.U or a string containing the name of an M-file which 
%   applies a linear operator to the columns of a given matrix.  In the latter
%   case, the second input argument must be N, the order of the problem.
%
%   If A is a struct with components A.L and A.U, such that 
%   L*U = (A - sigma*I), a shift-and-invert Lanczos iteration is performed
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
%     OPTIONS.Y    :  The lanczos vectors are reorthogonalized against
%                     the columns of the matrix OPTIONS.Y.
%
%   If both R0, Q_old and T_old are provided, they must contain 
%   a partial Lanczos tridiagonalization of A on the form
%
%        A Q_old = Q_old T_old + R0 .  
%
%   In this case the factorization is extended to dimension K x K by
%   continuing the Lanczos algorithm with R0 as starting vector.
%
%   On exit ANORM contains an approximation to ||A||_2. 
%     IERR = 0  :  K steps were performed succesfully.
%     IERR > 0  :  K steps were performed succesfully, but the algorithm
%                  switched to full reorthogonalization after IERR steps.
%     IERR < 0  :  Iteration was terminated after -IERR steps because an
%                  invariant subspace was found, and 3 deflation attempts 
%                  were unsuccessful.
%   On exit WORK(1) contains the number of reorthogonalizations performed, and
%   WORK(2) contains the number of inner products performed in the
%   reorthogonalizations.
%
%   See also LANEIG, REORTH, COMPUTE_INT

% References: 
% R.M. Larsen, Ph.D. Thesis, Aarhus University, 1998.
%
% G. H. Golub & C. F. Van Loan, "Matrix Computations",
% 3. Ed., Johns Hopkins, 1996.  Chapter 9.
%
% B. N. Parlett, ``The Symmetric Eigenvalue Problem'', 
% Prentice-Hall, Englewood Cliffs, NJ, 1980.
%
% H. D. Simon, ``The Lanczos algorithm with partial reorthogonalization'',
% Math. Comp. 42 (1984), no. 165, 115--142.

% Rasmus Munk Larsen, DAIMI, 1998


% Check input arguments.
if nargin<1, error('Not enough input arguments.');  end
if isnumeric(A) | isstruct(A)
  if isnumeric(A)
    [m n] = size(A);
    if m~=n | ~isequal(A,A') | ~isreal(A)
      error('A must be real symmetric')
    end  
  elseif isstruct(A)
    [m n] = size(A.L);
  end
    
  if nargin<7 | isempty(T_k), 
    anorm = []; est_anorm=1; 
  else
    anorm = T_k; est_anorm=0; 
  end
  if nargin<6,  Q_k=[]; T_k=[]; else,  T_k = Q_k; Q_k = options; end
  if nargin<4 | isempty(r),  options = []; else,  options = r;  end
  if nargin<3 | isempty(kmax),  
    r = rand(n,1)-0.5;
  else
    r = kmax;
  end
  if nargin<2 | isempty(nin);  kmax = max(10,n/10); else,  kmax = nin;  end   
else
  if nargin<2
    error('Not enough input arguments.');
  end
  % Check input functions and parse to create an internal object
  % if an explicit expression is given.
  [A, msg] = fcnchk(A);
  if ~isempty(msg)
    error(msg);
  end  
  n = nin;
  if nargin<8 | isempty(anorm), anorm = []; est_anorm=1; else est_anorm=0; end
  if nargin<7,  Q_k=[]; T_k=[]; end
  if nargin<5 | isempty(options),  options = [];          end
  if nargin<4 | isempty(r),  r = rand(n,1)-0.5;   end
  if nargin<3 | isempty(kmax);  kmax = max(10,n/10); end
end
 
% Set options.  
delta = sqrt(eps/kmax); % Desired level of orthogonality.
eta = eps^(3/4)/sqrt(kmax);     % Level of orth. after reorthogonalization.
cgs = 0;                % Flag for switching between iterated CGS and MGS.
elr = 1;                % Flag for switching extended local 
                        % reorthogonalization on and off.
deflate = 0;              % Flag for deflation against OPTIONS.Y
			
% Parse options struct
if ~isempty(options) & isstruct(options)
  c = fieldnames(options);
  for i=1:length(c)
    if strmatch(c(i),'delta'), delta = getfield(options,'delta');  end
    if strmatch(c(i),'eta'), eta = getfield(options,'eta'); end
    if strmatch(c(i),'cgs'), cgs = getfield(options,'cgs'); end
    if strmatch(c(i),'elr'), elr = getfield(options,'elr'); end
    if strmatch(c(i),'Y'), deflate = ~isempty(options.Y);  end
  end
end

np = 0;  nr = 0; ierr=0;

% Rule-of-thumb estimate on the size of round-off terms:
eps1 = sqrt(n)*eps/2; % Notice that {\bf u} == eps/2.
gamma = 1/sqrt(2);

% Prepare Lanczos iteration
if isempty(Q_k) % New Lanczos tridiagonalization.
  % Allocate space 
  alpha = zeros(kmax+1,1);  beta = zeros(kmax+1,1);
  Q_k = zeros(n,kmax);
  q = zeros(n,1); beta(1)=norm(r);
  omega = zeros(kmax,1); omega_max = omega;  omega_old = omega;
  omega(1) = 0;   force_reorth= 0;  
  j0 = 1;
else            % Extending existing Lanczos tridiagonalization.
  j = size(Q_k,2); % Size of existing factorization
  % Allocate space
  Q_k = [Q_k zeros(n,kmax-j)]; 
  alpha = zeros(kmax+1,1);  beta = zeros(kmax+1,1);
  alpha(1:j) = diag(T_k);  
  if j>1
    beta(2:j) = diag(T_k,-1);
  end
  q = Q_k(:,j);
  % Reorthogonalize r.
  beta(j+1) = norm(r);
  if j<kmax & beta(j+1)*delta < anorm*eps1,
    fro = 1;
  end
  if isfinite(delta)
    int = 1:j;
    [r,beta(j+1),rr] = reorth(Q_k,r,beta(j+1),int,gamma,cgs);
    np = rr*j;    nr = 1;   force_reorth = 1;  
  else
     force_reorth = 0;  
  end
  % Compute Gerscgorin bound on ||T_k||_2 as SQRT(||T_k'*T_k||_1)
  if est_anorm
    anorm = sqrt(norm(T_k'*T_k,1));
  end
  omega = eps1*ones(kmax,1); omega_max = omega;  omega_old = omega;
  j0 = j+1;
end

if delta==0
  fro = 1; % The user has requested full reorthogonalization.
else
  fro = 0;
end

for j=j0:kmax,  
  % Lanczos Step:
  q_old = q;
  if beta(j)==0
    q = r;
  else
    q = r / beta(j);
  end
  Q_k(:,j) = q;
  if isnumeric(A)
    u = A*q;
  elseif isstruct(A)
    u = A.U \ ( A.L \ q);
  else
    u = feval(A,q);
  end
  r = u - beta(j)*q_old;
  alpha(j) = q'*r;
  r = r - alpha(j)*q;
  

  % Extended local reorthogonalization:
  beta(j+1) = sqrt(r'*r); % Quick and dirty estimate.
  if beta(j+1)<gamma*beta(j) & elr 
    if  j==1
      t1=0;
      for i=1:2
	t = q'*r;    
	r = r-q*t;
	t1 = t1+t;
      end
      alpha(j) = alpha(j) + t1;
    elseif j>1
      t1 = q_old'*r;
      t2 = q'*r;
      r = r  - (q_old*t1 + q*t2); % Add small terms together first to
      if beta(j)~=0               % reduce risk of cancellation.
	beta(j) = beta(j) + t1;
      end
      alpha(j) = alpha(j) + t2;
    end        
    beta(j+1) = sqrt(r'*r); % Quick and dirty estimate.
  end

  % Update Gersgorin estimate of ||T_k|| if required
%  if est_anorm & beta(j+1)~=0
%    T_k = spdiags([[beta(2:j);0] alpha(1:j) beta(1:j)],-1:1,j,j);
%    anorm = sqrt(norm(T_k'*T_k,1))
%  end
  if  est_anorm & beta(j+1)~=0
    anorm = update_gbound(anorm,alpha,beta,j);
  end

  % Update omega-recurrence
  if j>1 & ~fro & beta(j+1)~=0
    [omega,omega_old] = update_omega(omega,omega_old,j,alpha,beta,...
	eps1,anorm);
    omega_max(j) = max(abs(omega));
  end

  % Reorthogonalize if required
  if j>1 & (fro  | force_reorth | omega_max(j)>delta) & beta(j+1)~=0
    if fro
      int = 1:j;
    else
      if force_reorth == 0
	force_reorth= 1; % Do forced reorth to avoid spill-over from q_{j-1}.
	int = compute_int(omega,j,delta,eta,0,0,0);
      else
	force_reorth= 0; 
      end
    end
    [r,beta(j+1),rr] = reorth(Q_k,r,beta(j+1),int,gamma,cgs);
    omega(int) = eps1;
    np = np + rr*length(int(:));    nr = nr + 1;
  else
    beta(j+1) = norm(r); % compute norm accurately.
  end

  if deflate    
    [r,beta(j+1),rr] = reorth(options.Y,r,beta(j+1),1:size(options.Y,2), ...
			      gamma,cgs);
  end
  
  if  j<kmax & beta(j+1) < n*anorm*eps  , 
    % If beta is "small" we deflate by setting the off-diagonals of T_k
    % to 0 and attempt to restart with a basis for a new 
    % invariant subspace by replacing r with a random starting vector:
    beta(j+1) = 0;
    bailout = 1;
    for attempt=1:3    
      r = rand(n,1)-0.5;  
      if isnumeric(A)
	r = A*r;
      elseif isstruct(A)
	r = A.U \ ( A.L \ r);
      else
	r = feval(A,r);
      end      
      nrm=sqrt(r'*r); % not necessary to compute the norm accurately here.
      int = 1:j;
      [r,nrmnew,rr] = reorth(Q_k,r,nrm,int,gamma,cgs);
      omega(int) = eps1;
      np = np + rr*length(int(:));    nr = nr + 1;
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
      r=r/nrmnew; % Continue with new normalized r as starting vector.
      force_reorth = 1;
      if delta>0
	fro = 0;    % Turn off full reorthogonalization.
      end
    end    
  elseif j<kmax & ~fro & beta(j+1)*delta < anorm*eps1,
    % If anorm*eps1/beta(j+1) > delta then  omega(j+1) will 
    % immediately exceed delta, and thus forcing a reorth. to occur at the
    % next step. The components of omega will mainly be determined
    % by the initial value and not the recurrence, and therefore we 
    % cannot tell reliably which components exceed eta => we might 
    % as well switch to full reorthogonalization to avoid trouble.
    % The user is probably trying to determine pathologically
    % small ( < sqrt(eps)*||A||_2 ) eigenvalues. 
    %    warning(['Semiorthogonality cannot be maintained at iteration ', ...
    %	  num2str(j),'. The matrix is probably ill-conditioned.', ...
    %	  ' Switching to full reorthogonalization.'])
    fro = 1;
    ierr = j;
  end
end

% Set up tridiagonal T_k in sparse matrix data structure.
T_k = spdiags([[beta(2:j);0] alpha(1:j) beta(1:j)],-1:1,j,j);
if nargout<2
  Q_k = T_k;
elseif j~=size(Q_k,2)
  Q_k = Q_k(:,1:j);
end
work = [nr np];


function [omega,omega_old] = update_omega(omega, omega_old, j, ...
    alpha,beta,eps1,anorm)
% UPDATE_OMEGA:  Update Simon's omega_recurrence for the Lanczos vectors.
%
% [omega,omega_old] = update_omega(omega, omega_old,j,eps1,alpha,beta,anorm)
% 

% Rasmus Munk Larsen, DAIMI, 1998.

% Estimate of contribution to roundoff errors from A*v 
%   fl(A*v) = A*v + f, 
% where ||f|| \approx eps1*||A||.
% For a full matrix A, a rule-of-thumb estimate is eps1 = sqrt(n)*eps.
T = eps1*anorm;
binv = 1/beta(j+1);

omega_old = omega;
% Update omega(1) using omega(0)==0.
omega_old(1)= beta(2)*omega(2)+ (alpha(1)-alpha(j))*omega(1) -  ...
    beta(j)*omega_old(1);
omega_old(1) = binv*(omega_old(1) + sign(omega_old(1))*T);
% Update remaining components.
k=2:j-2;
omega_old(k) = beta(k+1).*omega(k+1) + (alpha(k)-alpha(j)).*omega(k) ...
     + beta(k).*omega(k-1) - beta(j)*omega_old(k);
omega_old(k) = binv*(omega_old(k) + sign(omega_old(k))*T);       
omega_old(j-1) = binv*T;
% Swap omega and omega_old.
temp = omega;
omega = omega_old;
omega_old = omega;
omega(j) =  eps1;


function anorm = update_gbound(anorm,alpha,beta,j)
%UPDATE_GBOUND   Update Gerscgorin estimate of 2-norm 
%  ANORM = UPDATE_GBOUND(ANORM,ALPHA,BETA,J) updates the Gerscgorin bound
%  for the tridiagonal in the Lanczos process after the J'th step.
%  Applies Gerscgorins circles to T_K'*T_k instead of T_k itself
%  since this gives a tighter bound.

if j==1 % Apply Gerscgorin circles to T_k'*T_k to estimate || A ||_2
  i=j; 
  % scale to avoid overflow
  scale = max(abs(alpha(i)),abs(beta(i+1)));
  alpha(i) = alpha(i)/scale;
  beta(i+1) = beta(i+1)/scale;
  anorm = 1.01*scale*sqrt(alpha(i)^2+beta(i+1)^2 + abs(alpha(i)*beta(i+1)));
elseif j==2
  i=1;
  % scale to avoid overflow
  scale = max(max(abs(alpha(1:2)),max(abs(beta(2:3)))));
  alpha(1:2) = alpha(1:2)/scale;
  beta(2:3) = beta(2:3)/scale;
  
  anorm = max(anorm, scale*sqrt(alpha(i)^2+beta(i+1)^2 + ...
      abs(alpha(i)*beta(i+1) + alpha(i+1)*beta(i+1)) + ...
      abs(beta(i+1)*beta(i+2))));
  i=2;
  anorm = max(anorm,scale*sqrt(abs(beta(i)*alpha(i-1) + alpha(i)*beta(i)) + ...
      beta(i)^2+alpha(i)^2+beta(i+1)^2 +  ...
      abs(alpha(i)*beta(i+1))) );
elseif j==3
  % scale to avoid overflow
  scale = max(max(abs(alpha(1:3)),max(abs(beta(2:4)))));
  alpha(1:3) = alpha(1:3)/scale;
  beta(2:4) = beta(2:4)/scale;
  i=2;
  anorm = max(anorm,scale*sqrt(abs(beta(i)*alpha(i-1) + alpha(i)*beta(i)) + ...
      beta(i)^2+alpha(i)^2+beta(i+1)^2 +  ...
      abs(alpha(i)*beta(i+1) + alpha(i+1)*beta(i+1)) + ...
      abs(beta(i+1)*beta(i+2))) );
  i=3;
  anorm = max(anorm,scale*sqrt(abs(beta(i)*beta(i-1)) + ...
      abs(beta(i)*alpha(i-1) + alpha(i)*beta(i)) + ...
      beta(i)^2+alpha(i)^2+beta(i+1)^2 +  ...
      abs(alpha(i)*beta(i+1))) );
else
  % scale to avoid overflow
  %  scale = max(max(abs(alpha(j-2:j)),max(abs(beta(j-2:j+1)))));
  %  alpha(j-2:j) = alpha(j-2:j)/scale;
  %  beta(j-2:j+1) = beta(j-2:j+1)/scale;
  
  % Avoid scaling, which is slow. At j>3 the estimate is usually quite good
  % so just make sure that anorm is not made infinite by overflow.
  i = j-1;
  anorm1 = sqrt(abs(beta(i)*beta(i-1)) + ...
      abs(beta(i)*alpha(i-1) + alpha(i)*beta(i)) + ...
      beta(i)^2+alpha(i)^2+beta(i+1)^2 +  ...
      abs(alpha(i)*beta(i+1) + alpha(i+1)*beta(i+1)) + ...
      abs(beta(i+1)*beta(i+2)));
  if isfinite(anorm1)
    anorm = max(anorm,anorm1);
  end
  i = j;
  anorm1 = sqrt(abs(beta(i)*beta(i-1)) + ...
      abs(beta(i)*alpha(i-1) + alpha(i)*beta(i)) + ...
      beta(i)^2+alpha(i)^2+beta(i+1)^2 +  ...
      abs(alpha(i)*beta(i+1)));
  if isfinite(anorm1)
    anorm = max(anorm,anorm1);
  end
end
