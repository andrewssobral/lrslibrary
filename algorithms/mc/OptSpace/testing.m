clear;

% Sample Code to test the OptSpace program :

fprintf(1,'Generating Matrix...\n');
n = 1000;
m = 1000;
r = 3;
tol = 1e-8 ;

eps = 10*r*log10(n); 
% Generate the matrix with the given parameters
U = randn(n,r);
V = randn(m,r);
Sig = eye(r) ;
M0 = U*Sig*V';

% Select the entries independently with prob. eps/sqrt(mn) 
E = 1 - ceil( rand(n,m) - eps/sqrt(m*n)  ) ;
M_E = sparse(M0.*E) ;

fprintf(1,'Testing without noise...\n');
% Call the OptSpace function
[X S Y dist] = OptSpace(sparse(M_E),[],[],tol);

% Compute the Frobenius norm
err_nonoise = norm(X*S*Y' - M0,'fro')/sqrt(m*n);


% Add noise
sig = .01 ;
M = M0 + sig*randn(n,m);
M_E = sparse(M.*E) ;

fprintf(1,'Testing with noise...\n');
% Call the OptSpace function
[X S Y dist] = OptSpace(sparse(M_E),[],20,tol);

% Compute the Frobenius norm
err_noise = norm(X*S*Y' - M0,'fro')/sqrt(m*n);

fprintf(1,'RMSE (without noise)           : %e\n',err_nonoise);
fprintf(1,'RMSE (with noise var %f) : %e\n',sig,err_noise);

%%
numr = size(M,1);
numc = size(M,2);
I = randi([0 1],numr,numc);
MI = M.*I;
show_2dvideo(MI,m,n);
tol = 1e-8;
[X,S,Y,dist] = OptSpace(sparse(MI),[],20,tol);
L = X*S*Y';
show_2dvideo(L,m,n);

