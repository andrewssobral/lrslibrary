% Script for testing singular values computed using LANSVD, LANEIG and
% EIGS.

% Rasmus Munk Larsen, DAIMI, 1998.
global MxV MU MUTRUE NU NUTRUE
MU=[];
MUTRUE=[];
NU=[];
NUTRUE=[];
rand('state',0);
setup = 1; % set setup=1 to set up testproblem.
k=10;      % Number of singular values to compute
HB_dir = 'Harwell-Boeing';

% TESTPROBLEMS. Uncomment the relevant line to select a problem.
%===============================================================
%  NAME                 DIMENSIONS
%===============================================================

% Discrete 1. order derivative matrix:
%problem = 'derivative';%   500 x 501 

% Discrete 2. order derivative matrix. (symmetric)
%problem = 'laplace';   %   500 x 500

% Test problem from helioseismology.
%problem='helio212a';   %   212 x 100
%problem='helio212b';   %   212 x 100

% Harwell-Boeing testmatrices:
%problem='abb313';      %   313 x   175
%problem='west0479';    %   479 x   479
%problem='rbs480a';     %   480 x   480
%problem='illc1850';    %  1850 x   712
%problem='qh1484';      %  1484 x  1484
problem='mhd4800a';    %  4800 x  4800
%problem='cry10000';    % 10000 x 10000
%problem='fidap011';    % 16614 x 16614
%problem='af23560';     % 23560 x 23560
%problem='bcsstk32';    % 44609 x 44609 
%problem='s3dkt3m2';    % 90449 x 90449

global A C
format compact

if setup==1,
  A=[]; C=[]; Sigma=[];
  switch problem
    case 'derivative'
      m = 500; n=m+1;
      e = ones(n,1);
      A = spdiags([e -e], 0:1, n-1, n); 
      Sigma = -sort(-2*cos(pi/2*[m:-1:1]'./(m+1))); 
      Sigma1 = svd(full(A));
    case 'laplace'
      m = 100; n=m;
      e = ones(n,1);
      A = spdiags([-e 2*e -e], -1:1, n, n); 
      Sigma = -sort(-4*cos(pi/2*[n:-1:1]'./(n+1)).^2); 
      Sigma1 = svd(full(A));
    case 'helio212a'
      load 'helio.mat'
      [m n] = size(A);
      [U,Sigma,V] = svd(full(A));
      x = rand(m,1); x = x/norm(x);
      [dummy,P1] = sort(x);
      y = rand(n,1);    y = y/norm(y);
      [dummy,P2] = sort(y);
      A = sparse(Sigma);
      A = A(P1,P2); % Permute rows and columns randomly.
      Sigma = diag(Sigma);
      %      shift = eps*min(Sigma);      
      shift = 0;
      % Add a small perturbation to try hide that A is just
      % a permuted diagonal and see if we can fool SVD into doing full
      % reduction to bidiagonal form.
      if shift~=0
	A = A + shift*(x*y');
	Sigma = Sigma+shift;
      end      
      Sigma1 = svd(full(A));
    case 'helio212b'
      load 'helio.mat'
      [U,Sigma,V] = svd(full(A));
      A = U * Sigma *V';
      Sigma = diag(Sigma);
      Sigma1 = svd(full(A));
    case 'diag'      
      m = 200; n=100;
%      Sigma = logspace(0,-0.1,n);
      Sigma = [n:-1:1]';
      Sigma(1) = 1000;
      A = [spdiags(Sigma(:),0,n,n); sparse(m-n,n)];      
    otherwise
      cur = pwd;
      cd(HB_dir)
      file = [problem,'.mtx'];
      [A,m,n,nz] = mmread(file);      
      cd(cur)
      A = sparse(A);
      if min(m,n)^2*max(m,n) < 300^3
	Sigma = svd(full(A));
      end
  end
  [m n] = size(A);
  fprintf('nnz(A) = %i\n',nnz(A));
  C = sparse([[sparse(m,m),A];[A',sparse(n,n)]]);
end

%%% Print information about matrix %%%
clc;
disp('************** PROPACK TEST ****************')
fprintf('The testmatrix %s is %i-by-%i and has %i non-zero elements.\n', ...
    problem,m,n,nnz(A));
if issparse(A) 
  fprintf('Here is what the sparsity pattern of the matrix looks like... \n(press a key to continue)\n')
  spy(A)
  title(problem)  
  pause
hold off
end
fprintf('Calculating %i singular triplets of the matrix.\n',k)


if max(m,n) < 300
%%%%%%%%%%%%%%%%%%%% SVD %%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('-------------------------------------------------------------------\n')
fprintf('\n')
fprintf('METHOD =  SVD\n')
fprintf('\n')


flops(0);tic;
[U,S,V] = svd(full(A),0);
Sigma=diag(S);
fprintf('Elapsed time    = %f\n',toc);
fprintf('Number of flops = %e\n',flops);
fprintf('Residual norm   = %e\n',norm(A*V - U*S,'fro')/norm(A,'fro'));
fprintf('Orthogonality U = %e\n',norm(eye(n) - U'*U,'fro'))
fprintf('Orthogonality V = %e\n',norm(eye(n) - V'*V,'fro'))
end


%%%%%%%%%%%%%%%%%%%% LANSVD %%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('-------------------------------------------------------------------\n')
fprintf('\n')
fprintf('METHOD =  LANSVD\n')
fprintf('\n')

options = [];
v0 = rand(m,1)-0.5;
options.p0 = v0;

MxV = 0;flops(0);tic;
[U1,S1,V1] = lansvd('Afunc','Atransfunc',m,n,k,'L',options);
fprintf('Number MxV      = %f\n',MxV);
fprintf('Elapsed time    = %f\n',toc);
fprintf('Number of flops = %e\n',flops);
fprintf('Residual norm   = %e\n',norm(A*V1 - U1*S1,'fro')/norm(A,'fro'));
fprintf('Orthogonality U = %e\n',norm(eye(k) - U1'*U1,'fro'))
fprintf('Orthogonality V = %e\n',norm(eye(k) - V1'*V1,'fro'))

if ~isempty(Sigma)
  % Compare with "true" singular values
  S1 = diag(S1);
  l= length(S1);
  E1 = abs((S1 - Sigma(1:l))./Sigma(1:l));
  fprintf('Max rel.  error = %g\n',max(E1))
  fprintf('Mean rel. error = %g\n',mean(E1))
  if exist('Sigma1') & ~isempty(Sigma1)
    semilogy(abs((Sigma1(1:l) - Sigma(1:l))./Sigma(1:l)),'+');    
    hold on
  end
  semilogy(E1,'o');  
  hold off
  ylabel('Relative error  |\theta_i - \sigma_i| / \sigma_i')
  xlabel('i')
end

%%%%%%%%%%%%%%%%%%%% LANEIG(A'*A) %%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('-------------------------------------------------------------------\n')
fprintf('\n')
fprintf('METHOD =  LANEIG(A''*A)\n')
fprintf('\n')
options=[];
options.v0 = A'*v0;


flops(0);MxV=0; tic;
%S = laneig(C,k,'AL',options);
[V3,S3,b3] = laneig('AtAfunc',n,k,'AL',options);
U3 = A*V3;
for i=1:k
  U3(:,i) = U3(:,i)/norm(U3(:,i));
end

fprintf('Number MxV      = %f\n',MxV);
fprintf('Elapsed time    = %f\n',toc);
fprintf('Number of flops = %e\n',flops);
fprintf('Residual norm   = %e\n',norm(A'*(A*V3) - V3*S3,'fro')/norm(A'*A,'fro'));
fprintf('Orthogonality U = %e\n',norm(eye(k) - U3'*U3,'fro'))
fprintf('Orthogonality V = %e\n',norm(eye(k) - V3'*V3,'fro'))
S3 = diag(sqrt(diag(S3)));

if ~isempty(Sigma)
  S3 = diag(S3);
  l= length(S3);
  % Compare with "true" singular values
  E3 = abs((S3 - Sigma(1:l))./Sigma(1:l));
 fprintf('Max rel.  error = %g\n',max(E3))
 fprintf('Mean rel. error = %g\n',mean(E3))
  hold on
  semilogy(E3,'d');  
  hold off
end

%%%%%%%%%%%%%%%%%%%% LANEIG(C) %%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('-------------------------------------------------------------------\n')
fprintf('\n')
fprintf('METHOD =  LANEIG(C)\n')
fprintf('\n')
options.v0 = [v0;zeros(n,1)];

tic; flops(0);MxV=0;
%S = laneig(C,k,'AL',options);
[V2,S2] = laneig('Cfunc',m+n,k,'AL',options);
U2= sqrt(2)*V2(1:m,:);
V2= sqrt(2)*V2(m+1:end,:);
fprintf('Number MxV      = %f\n',MxV);
fprintf('Elapsed time    = %f\n',toc);
fprintf('Number of flops = %e\n',flops);
fprintf('Residual norm   = %e\n',norm(A*V2 - U2*S2,'fro')/norm(A,'fro'));
fprintf('Orthogonality U = %e\n',norm(eye(k) - U2'*U2,'fro'))
fprintf('Orthogonality V = %e\n',norm(eye(k) - V2'*V2,'fro'))


if ~isempty(Sigma)
  S2 = diag(S2);
  l= length(S2);
  % Compare with "true" singular values
  E2 = abs((S2 - Sigma(1:l))./Sigma(1:l));
  fprintf('Max rel.  error = %g\n',max(E2))
  fprintf('Mean rel. error = %g\n',mean(E2))
  hold on
  semilogy(E2,'x');  
  hold off
  legend('SVD','BPRO','PRO(A''*A)','PRO(C)',2)
end


%%%%%%%%%%%%%%%%%%%% EIGS %%%%%%%%%%%%%%%%%%%%%%%%
% Try ARPACK based eigs routine:
fprintf('\n')
fprintf('-------------------------------------------------------------------\n')
fprintf('METHOD =  EIGS(C)\n')
fprintf('\n')

options =[];
options.disp=0;
options.issym = 1;
options.v0 = [v0;zeros(n,1)];

tic; flops(0);MxV=0;
[V4,S4] = eigs('Cfunc',m+n,k,'LR',options);
U4= sqrt(2)*V4(1:m,:);
V4= sqrt(2)*V4(m+1:end,:);
fprintf('Number MxV      = %f\n',MxV);
fprintf('Elapsed time    = %f\n',toc);
fprintf('Number of flops = %e\n',flops);
fprintf('Residual norm   = %e\n',norm(A*V4 - U4*S4,'fro')/norm(A,'fro'));
fprintf('Orthogonality U = %e\n',norm(eye(k) - U4'*U4,'fro'))
fprintf('Orthogonality V = %e\n',norm(eye(k) - V4'*V4,'fro'))

if ~isempty(Sigma)
  % Compare with "true" singular values
  S4 = diag(S4);
  l= length(S4);
  E4 = abs((S4 - Sigma(1:l))./Sigma(1:l));
 fprintf('Max rel.  error = %g\n',max(E4))
 fprintf('Mean rel. error = %g\n',mean(E4))
  hold on
  semilogy(E4,'*');
  hold off
  legend('SVD','BPRO','PRO(A''*A)','PRO(C)','EIGS',2)
end



%%%%%%%%%%%%%%%%%%%% EIGS %%%%%%%%%%%%%%%%%%%%%%%%
% Try ARPACK based eigs routine:
fprintf('-------------------------------------------------------------------\n')
fprintf('\n')
fprintf('METHOD =  EIGS(A''*A)\n')
fprintf('\n')

options =[];
options.disp=0;
options.issym = 1;
options.v0 = A'*v0;

tic; flops(0); MxV=0;
[V5,S5] = eigs('AtAfunc',n,k,'LR',options);
U5 = A*V5;
for i=1:k
  U5(:,i) = U5(:,i)/norm(U5(:,i));
end
fprintf('Number MxV      = %f\n',MxV);
fprintf('Elapsed time    = %f\n',toc);
fprintf('Number of flops = %e\n',flops);
fprintf('Residual norm   = %e\n',norm(A'*(A*V5) - V5*S5,'fro')/norm(A'*A,'fro'));
fprintf('Orthogonality U = %e\n',norm(eye(k) - U5'*U5,'fro'))
fprintf('Orthogonality V = %e\n',norm(eye(k) - V5'*V5,'fro'))
S5 = diag(sqrt(diag(S5)));

if ~isempty(Sigma)
  % Compare with "true" singular values
  S5 = diag(S5);
  l= length(S5);
  E5 = abs((S5 - Sigma(1:l))./Sigma(1:l));
  fprintf('Max rel.  error = %e\n',max(E5))
  fprintf('Mean rel. error = %e\n',mean(E5))
  hold on
  semilogy(E5,'^');
  hold off
  title('Accuracy of \sigma_i compared with output from the svd command.')
  legend('BPRO','PRO(A''*A)','PRO(C)','EIGS(C)','EIGS(A''*A)',2)
end
