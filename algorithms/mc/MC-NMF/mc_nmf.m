function [X,Y,Out] = mcnf(A,Omega,k,m,n,opts)
%%%
% solver for matrix completion with nonnegative factors
%
%        min ||P_omega(XY - M)||_F^2, 
%        st.    X,Y >= 0
%
% Output:
%           X --- m x k matrix
%           Y --- k x n matrix
%          Out --- output information
% Input:
%           A --- given partial matrix M(Omega)
%           k --- given estimate rank
%           opts --- option structure
%
% Copyright(c) 2011 Yangyang Xu
%
%
% Based on the code of lmafit_nmf.m by Yin Zhang
%        and lmafit_mc_adp.m by Yin Zhang and Zaiwen Wen
% In the function, two other functions used in lmafit_mc_adp.m 
% are called for sparse case

if min(A) < 0
    error('Matrix entries must be nonnegative');
end

L = length(A);
Zfull = (L/(m*n) > 0.2 ) || k > .02*min(m,n) || m*n < 5e5;

%set parameters;
tol = 5e-6; maxit = 500; iprint = 1; Y0 = rand(k,n); gamma = 1.618;
        
if isfield(opts,'tol');       tol = opts.tol;    end
if isfield(opts,'maxit');   maxit = opts.maxit;  end
if isfield(opts,'print');  iprint = opts.print;  end
if isfield(opts,'Y0');         Y0 = opts.Y0;     end
if isfield(opts,'Zfull');   Zfull = opts.Zfull;  end

linopts.SYM = true;    linopts.POSDEF = true;
relres = 1;  nstall = 0;  I = speye(k);       

%scale the problem, different value of Mnrm can be used
Mnrm = 2.5e+5;
scal = Mnrm/norm(A);
A = scal*A; 
r1 = max(m,n) / min(m,n);
r2 = k / min(m,n);

%set penalty parameters
beta1 = Mnrm* (r1 / r2) * 2e-4;
beta2 = (n/m)*beta1;

if iprint == 2; fprintf('Initial beta: %6.3e%6.3e\n',beta1, beta2); end
if iprint == 1; fprintf('Iteration:     '); end

% initialize
X = zeros(m,k); Y = Y0;
Y = Y / norm(Y,'fro') * sqrt(Mnrm);
U = zeros(m,k); Lam1 = zeros(m,k);
V = zeros(k,n); Lam2 = zeros(k,n);

if Zfull %Z is full
    Z = zeros(m,n);  Z(Omega) = A;   
else %Z = S + XY, initialize the storage of S
    A(A==0) = eps;
    if isnumeric(Omega);       [Ik,Jk] = ind2sub([m n],Omega);
    elseif isstruct(Omega)     Ik = Omega.Ik; Jk = Omega.Jk;     end
    %make sure the order of Ik, Jk and data are correctly as in a sparse matrix
    S = sparse(Ik, Jk, A, m, n); [Ik, Jk, A] = find(S);   A = A';
end


for iter = 1:maxit
      
    % updating variables X, U and Lam1
    gX = beta1*U - Lam1;
    if Zfull     Xt = (gX + Z*Y')';
    else     Yt = Y'; Xt = (gX + S*Yt + X*(Y*Yt))';
    end
    B = Y*Y' + beta1*I;
    Xt = linsolve(B, Xt, linopts);
    X = Xt'; 
    if Zfull    XtZ = Xt*Z;
    else      XtZ = Xt*S+(Xt*X)*Y;
    end
    U = max(0, X + Lam1/beta1);
    Lam1 = Lam1 + gamma*beta1*(X-U);
    
    % updating variables Y, V and Lam2
    gY = beta2*V - Lam2; 
    B = Xt*X + beta2*I;
    Y = linsolve(B, gY+XtZ, linopts);
    V = max(0, Y + Lam2/beta2);    
    Lam2 = Lam2 + gamma*beta2*(Y-V);
    
    % updating variable Z
    if Zfull   Z = X*Y;  Res = A-Z(Omega);  Z(Omega) = A;
    else
        Res = A-partXY(Xt, Y, Ik, Jk, L);
        updateSval(S, Res, L);
    end
    
    relres0 = relres;
    if Zfull    relres=norm(Res)/Mnrm;
    else    relres = norm(Res)/Mnrm;
    end

    % printout
    if iprint == 1; fprintf('\b\b\b\b\b%5i',iter); end
    if iprint == 2;        
        fprintf('iter %5i:  relres %6.3e \n',...
            iter,relres);
    end
    
    % check stopping
    crit1 = abs(relres0-relres) < tol*max(1,relres0);
    crit2 = relres < tol;
    
    if crit1; nstall = nstall + 1; else nstall = 0; end; 
    %make crit1 satisfy in three successive times
    
    if (crit1 || crit2)  && nstall > 3
        if iprint == 2;
            fprintf('crits = [%i %i]\n',...
                crit1,crit2);
        end
        break; 
    end
    
end %iter

X = max(0,X)/sqrt(scal); 
Y = max(0,Y)/sqrt(scal); 

if iprint == 1; fprintf('\n'); end

Out.iter = iter;
Out.Lam1 = Lam1;
Out.Lam2 = Lam2;
Out.beta = [beta1 beta2];

end%main