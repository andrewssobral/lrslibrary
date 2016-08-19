function [V,D,bnd,j,work] = laneig(A,nin,k,sigma,options)

%LANEIG  Compute a few eigenvalues and eigenvectors.
%   LANEIG solves the eigenvalue problem A*v=lambda*v, when A is 
%   real and symmetric using the Lanczos algorithm with partial 
%   reorthogonalization (PRO). 
%
%   [V,D] = LANEIG(A) 
%   [V,D] = LANEIG('Afun',N) 
%
%   The first input argument is either a real symmetric matrix, or a 
%   string containing the name of an M-file which applies a linear 
%   operator to the columns of a given matrix.  In the latter case,
%   the second input argument must be N, the order of the problem.
%
%   The full calling sequence is
%
%   [V,D,ERR] = LANEIG(A,K,SIGMA,OPTIONS)
%   [V,D,ERR] = LANEIG('Afun',N,K,SIGMA,OPTIONS)
%
%   On exit ERR contains the computed error bounds.  K is the number of
%   eigenvalues desired and SIGMA is numerical shift or a two letter string
%   which specifies which part of the spectrum should be computed:
%
%   SIGMA            Specified eigenvalues
%
%   'AL'            Algebraically Largest 
%   'AS'            Algebraically Smallest
%   'LM'            Largest Magnitude   (default)
%   'SM'            Smallest Magnitude  (does not work when A is an m-file)
%   'BE'            Both Ends.  Computes k/2 eigenvalues
%                   from each end of the spectrum (one more
%                   from the high end if k is odd.) 
%
%   The OPTIONS structure specifies certain parameters in the algorithm.
%
%    Field name      Parameter                              Default
%   
%    OPTIONS.tol     Convergence tolerance                  16*eps
%    OPTIONS.lanmax  Dimension of the Lanczos basis.
%    OPTIONS.v0      Starting vector for the Lanczos        rand(n,1)-0.5
%                    iteration.
%    OPTIONS.delta   Level of orthogonality among the       sqrt(eps/K)
%                    Lanczos vectors.
%    OPTIONS.eta     Level of orthogonality after           10*eps^(3/4)
%                    reorthogonalization. 
%    OPTIONS.cgs     reorthogonalization method used        0
%                    '0' : iterated modified Gram-Schmidt 
%                    '1' : iterated classical Gram-Schmidt
%    OPTIONS.elr     If equal to 1 then extended local      1
%                    reorthogonalization is enforced. 
%
%   See also LANPRO, EIGS, EIG.

% References: 
% R.M. Larsen, Ph.D. Thesis, Aarhus University, 1998.
%
% B. N. Parlett, ``The Symmetric Eigenvalue Problem'', 
% Prentice-Hall, Englewood Cliffs, NJ, 1980.
%
% H. D. Simon, ``The Lanczos algorithm with partial reorthogonalization'',
% Math. Comp. 42 (1984), no. 165, 115--142.

% Rasmus Munk Larsen, DAIMI, 1998


%%%%%%%%%%%%%%%%%%%%% Parse and check input arguments. %%%%%%%%%%%%%%%%%%%%%%

if ~isstr(A)
  if nargin<1
    error('Not enough input arguments.');
  end
  [m n] = size(A);
  Aisfunc = 0;
  if m~=n | ~isequal(A,A') | ~isreal(A)
    error('A must be real symmetric')
  end  
  if nargin < 4 | isempty(sigma)
    options = [];
  else  
    options = sigma; 
  end
  if nargin < 3 | isempty(k), sigma = 'LM'; else, sigma = k; end
  if nargin < 2 | isempty(nin), k = min(n,5); else, k = nin; end
else
  if nargin<2
    error('Not enough input arguments.');
  end
  Aisfunc = 1;
  n = nin;
  if nargin < 5 | isempty(options)
    options.tol = 16*eps;
    options.lanmax = n;
    options.v0 = rand(n,1)-0.5;
  end
  if nargin < 4 | isempty(sigma), sigma = 'LM'; end
  if nargin < 3 | isempty(k), k = min(n,5);  end
end

if ~isnumeric(k) | real(abs(fix(k)))~=k | ~isnumeric(n) | real(abs(fix(n)))~=n
  error('Input arguments N and K must be positive integers.')
end

% Quick return for n<2  or k<1
if n < 1 | k<1
  if nargout < 2
    V = zeros(k,1);
  else
    V = eye(n,k);
    D = zeros(k,k);
    bnd =zeros(k,1);
  end
  return
end
if n == 1 
  if ~Aisfunc
    D = A;
    V = 1;
    bnd = 0;
  else
    D = feval(A,1);
    V = 1;
    dnb = 0;
  end
  if nargout<2
    V=D;
  end
  return
end

% A is the matrix of all zeros (not detectable if A is a string)
if ~Aisfunc 
  if nnz(A)==0
    if nargout < 2
      V = zeros(k,1);
    else
      V = eye(n,k);
      D = zeros(k,k);
      bnd =zeros(k,1);
    end
    return
  end
end

lanmax = n;
tol = 16*eps;
r = rand(n,1)-0.5;
part = sigma;
% Parse options struct
if ~isempty(options) & isstruct(options)
  c = fieldnames(options);
  for i=1:length(c)
    if strmatch(c(i),'v0'), r = getfield(options,'v0'); r=r(:); end
    if strmatch(c(i),'tol'), tol = getfield(options,'tol'); end
    if strmatch(c(i),'lanmax'), lanmax = getfield(options,'lanmax'); end
  end
end

% Protect against absurd arguments.
tol = max(tol,eps);
lanmax = min(lanmax,n);
if size(r,1)~=n
  error('v0 must be a vector of length n')
end

lanmax = min(lanmax,n);
if k>lanmax
  error('K must satisfy  K <= LANMAX <= N.');
end
ksave = k;

if strcmp(sigma,'SM') & ~isstr(A)
  sigma = 0;
end


% Prepare for shift-and-invert if sigma is numeric.
if  isnumeric(sigma)
  part = 'LM';
  if isstr(A) 
    error('Shift-and-invert works only when the matrix A is given explicitly.');
  else
    pmmd = symmmd(A);
    A = A(pmmd,pmmd);
    [S.L,S.U] = lu(A - sigma*speye(n));
    condU = condest(S.U);
    dsigma = n * full(max(max(abs(A)))) * eps;
    if sigma < 0
      sgnsig = -1;
    else
      sgnsig = 1;
    end
    sigitr = 1;
    while condU > 1/eps & ((dsigma <= 1 & sigitr <= 10) | ~isfinite(condU))
      disps1 = sprintf(['sigma = %10e is near an exact eigenvalue of A,\n' ...
			'so we cannot use the LU factorization of (A-sigma*I): ' ...
			' condest(U) = %10e.\n'],sigma,condU);
      if abs(sigma) < 1
	sigma = sigma + sgnsig * dsigma;
	disps2 = sprintf('We are trying sigma + %10e = %10e instead.\n', ...
			 sgnsig*dsigma,sigma);
      else
	sigma = sigma * (1 + dsigma);
	disps2 = sprintf('We are trying sigma * (1 + %10e) = %10e instead.\n', ...
			 dsigma,sigma);
      end
      %     if nargout < 3 & dispn ~= 0             
      disp([disps1 disps2])
      %     end   
      [S.L,S.U] = lu(A - sigma*speye(n));
      condU = condest(S.U);
      dsigma = 10 * dsigma;
      sigitr = sigitr + 1;
    end
  end
  A = S;
end


neig = 0; nrestart=-1;
if ~strcmp(part,'BE') 
  j = min(2*k+2,lanmax);
else
  j = min(k+1,lanmax);
end


%%%%%%%%%%%%%%%%%%%%% Here begins the computation  %%%%%%%%%%%%%%%%%%%%%%

V = []; T = []; anorm = []; work = zeros(1,2); rnorm=-1;




while neig < k 
  %%%%%%%%%%%%%%%%%%%%% Compute Lanczos tridiagonalization %%%%%%%%%%%%%%%%%
  j = min(lanmax,j+1-mod(j,2));
  % "Trick" to avoid unwanted zero eigenvalues when laneig is used for
  % SVD calculations. (Nothing to if lanmax is odd, though.)
  
  if  ~isstr(A)
    [V,T,r,anorm,ierr,w] = lanpro(A,j,r,options,V,T,anorm);
  else
    [V,T,r,anorm,ierr,w] = lanpro(A,n,j,r,options,V,T,anorm);
  end
  work= work + w;

  if ierr<0 % Invariant subspace of dimension -ierr found. 
    j = -ierr;
  end

  %%%%%%%%%%%%%%%%%% Compute eigenvalues and error bounds %%%%%%%%%%%%%%%%%%
  % Analyze T
  [D,top,bot,err] = tqlb([full(diag(T))],full([0;diag(T,1)]));
  %  if err>0
  %    printf(['TQLB failed. Eigenvalue no. %i did not converge in 30', ...
  %	  ' iterations'],err);
  %  end
  %  full(T)
  %  [P,D] = eig(full(T));
  %  D = diag(D);
  %  bot = P(end,:)';
  %  [P(1,:)' P(end,:)']
  [D,I] = sort(D);
  bot = bot(I);
  
  % Set simple error bounds
  rnorm = norm(r);
  bnd = rnorm*abs(bot);
  
  % Use Largest Ritz value to estimate ||A||_2. This might save some
  % reorth. in case of restart.
  anorm = max(abs(D));
  
  % Estimate gap structure and refine error bounds
  bnd = refinebounds(D,bnd,n*eps*anorm);

  %%%%%%%%%%%%%%%%%%% Check convergence criterion %%%%%%%%%%%%%%%%%%%%
  % Reorder eigenvalues according to SIGMA
  switch part
   case 'AS'
    IPART = 1:j;
   case 'AL' 
    IPART = j:-1:1;
   case 'LM'
    [dummy,IPART] = sort(-abs(D));
   case 'BE'
    if j<k
      IPART=1:j;
    else
      mid = floor(k/2);
      par = rem(k,1);
      IPART = [1:mid,(j-mid-par):j]';
    end    
   otherwise
    error(['Illegal value for SIGMA: ',part]);
  end
  D = D(IPART);  bnd = bnd(IPART);
  if isnumeric(sigma)
    D = sigma + 1./D;
  end
  
  % Check if enough have converged.
  neig = 0;
  for i=1:min(j,k)
    if bnd(i) <= tol*abs(D(i))
      neig = neig + 1;
    end
  end
  
  %%%%%%%%%%% Check whether to stop or to extend the Krylov basis? %%%%%%%%%%
  if ierr<0 % Invariant subspace found
    if j<k
      warning(['Invariant subspace of dimension ',num2str(j-1),' found.'])
    end
    break;
  end
  if j>=lanmax % Maximal dimension of Krylov subspace reached => Bail out!
    if neig<ksave
      warning(['Maximum dimension of Krylov subspace exceeded prior',...
	       ' to convergence.']);
    end
    break;
  end
  
  % Increase dimension of Krylov subspace and try again.
  if neig>0
    %    j = j + ceil(min(20,max(2,((j-1)*(k-neig+1))/(2*(neig+1)))));
    j = j + min(100,max(2,0.5*(k-neig)*j/(neig+1)));
  elseif neig<k
    %    j = j + ceil(min(20,max(8,(k-neig)/2)));
    j = max(1.5*j,j+10);
  end
  j = min(j+1,lanmax);
  nrestart = nrestart + 1;
end



%%%%%%%%%%%%%%%% Lanczos converged (or failed). Prepare output %%%%%%%%%%%%%%%
k = min(ksave,j);

if nargout>1
  j = size(T,1);
  [Q,D] = eig(full(T)); D = diag(D);
  [D,I] = sort(D);
  % Compute and normalize Ritz vectors (overwrite V to save memory).
  V = V*Q(:,I(IPART(1:k)));
  for i=1:k
    nq = norm(V(:,i));
    if isfinite(nq) & nq~=0 & nq~=1
      V(:,i) = V(:,i)/nq;
    end
  end
  [D,I] = sort(D);
  D = D(IPART(1:k));
  if isnumeric(sigma)
    D = sigma + 1./D;
    V(pmmd,:) = V;
  end
end

% Pick out desired part of the spectrum
if length(D)~=k
  D = D(1:k);
  bnd = bnd(1:k);
end

if nargout<2
  V = D;
else
  D = diag(D);
end

