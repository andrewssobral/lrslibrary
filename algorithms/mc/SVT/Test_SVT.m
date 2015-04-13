%{
    Thanks for downloading SVT.

    SVT, and PROPACK (a software package by R.M.Larsen for sparse
        SVDs) contain some .mex files.  These files have already
        been compiled for most standard computer platforms, so
        you probably do not need to do anything special.  But if
        you need to compile them, go the "private" subdirectory
        and run "install_mex".

    Also in the "private" subdirectory is "test_PROPACK.m", which
        tests the PROPACK installation.  We have included only the
        part of PROPACK relevant to SVT (and made a few small changes).

    This files runs SVT to verify that it works.

    -Stephen Becker, March 2009  srbecker@caltech.edu

    Major changelog:
        5/13/09 should be able to handle complex numbers,
            though this hasn't been extensively tested.
            PROPACK, SVT, and smvp were modified to allow complex data.
        10/2/09 small change to a parameter after it was discovered
            that PROPACK wasn't returning accurate results on
            some platforms, causing lack of convergence.

%}
if ispc
    % reorth.f isn't compiled for Windows, but this shouldn't be an
    % issue, because the .m file version is pretty fast
    warning('off','PROPACK:NotUsingMex');
end
%% Setup a matrix
randn('state',2009);
rand('state',2009);

% if you're daring, try with complex numbers:
COMPLEX = false;

n1 = 150; n2 = 300; r = 10;
if COMPLEX
    M = (randn(n1,r)+1i*randn(n1,r))*(randn(r,n2)+1i*randn(r,n2))/2;
else
    M = randn(n1,r)*randn(r,n2);
end

df = r*(n1+n2-r);
oversampling = 5; 
m = min(5*df,round(.99*n1*n2) ); 
p  = m/(n1*n2);

Omega = randsample(n1*n2,m);  % this requires the stats toolbox
% a workaround, if you don't have the stats toolbox, is this:
%   Omega = randperm(n1*n2); Omega = Omega(1:m);

data = M(Omega);
% add in noise, if desired
sigma = 0;
% sigma = .05*std(data);
data = data + sigma*randn(size(data));

fprintf('Matrix completion: %d x %d matrix, rank %d, %.1f%% observations\n',...
    n1,n2,r,100*p);
fprintf('\toversampling degrees of freedom by %.1f; noise std is %.1e\n',...
    m/df, sigma );
if ~isreal(M), disp('Matrix is complex'); end
%% Set parameters and solve

tau = 5*sqrt(n1*n2); 
delta = 1.2/p;    
%{
 if n1 and n2 are very different, then
   tau should probably be bigger than 5*sqrt(n1*n2)

 increase tau to increase accuracy; decrease it for speed

 if the algorithm doesn't work well, try changing tau and delta
   i.e. if it diverges, try a smaller delta (e.g. delta < 2 is a 
   safe choice, but the algorithm may be slower than necessary).
%}
maxiter = 500; 
tol = 1e-4;
%% Approximate minimum nuclear norm solution by SVT algorithm
% Note: SVT, as called below, is setup for noiseless data 
%   (i.e. equality constraints).

fprintf('\nSolving by SVT...\n');
tic
[U,S,V,numiter] = SVT([n1 n2],Omega,data,tau,delta,maxiter,tol);
toc
    
X = U*S*V';
    
% Show results
fprintf('The recovered rank is %d\n',length(diag(S)) );
fprintf('The relative error on Omega is: %d\n', norm(data-X(Omega))/norm(data))
fprintf('The relative recovery error is: %d\n', norm(M-X,'fro')/norm(M,'fro'))
fprintf('The relative recovery in the spectral norm is: %d\n', norm(M-X)/norm(M))

%% Approximate minimum nuclear norm solution by FPC algorithm
% This version of FPC uses PROPACK for the multiplies
%   It is not optimized, and the parameters have not been tested
% The version in the FPC paper by Shiqian Ma, Donald Goldfarb and Lifeng Chen
%   uses an approximate SVD that will have different properties;
%   That code may be found at http://www.columbia.edu/~sm2756/FPCA.htm

mu_final = .01; tol = 1e-3;

fprintf('\nSolving by FPC...\n');
tic
[U,S,V,numiter] = FPC([n1 n2],Omega,data,mu_final,maxiter,tol);
toc
   
X = U*S*V';

% Show results
fprintf('The recovered rank is %d\n',length(diag(S)) );
fprintf('The relative error on Omega is: %d\n', norm(data-X(Omega))/norm(data))
fprintf('The relative recovery error is: %d\n', norm(M-X,'fro')/norm(M,'fro'))
fprintf('The relative recovery in the spectral norm is: %d\n', norm(M-X)/norm(M))
    
