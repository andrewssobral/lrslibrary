function [U,S,V,numiter]  = FPC(n,Omega,b,mu_final,maxiter,tol)
% [U,Sigma,V,numiter]  = FPC(n,Omega,b,mu_final,maxiter,gtol)
%
% Finds mininum  mu ||X||_* + 1/2 || A(X) - b ||_2^2
%
%   where A(X) is the projection of X onto the set of indices
%   in Omega.
%
% For efficiency, the algorithm uses continuation (i.e. a series of
%   mu, aka the "outer loop"), until mu = mu_final.
%
% maxiter controls maximum number of iterations per inner loop
%
% Outputs:
%   U,V and Sigma are singular vectors and singular values, where
%       X = U*Sigma*V'
%   numiter is the number of iterations over all inner and outer loops

% Reference:
%
%   "Fixed point and Bregman iterative methods for matrix
%       rank minimization."
%   Shiquian Ma, Donald Goldfarb, Lifeng Chen, October 2008
%   ftp://ftp.math.ucla.edu/pub/camreport/cam08-78.pdf

% code by Stephen Becker, srbecker@caltech.edu, March 2009

% May 2009: adding support for complex matrices


% -- some parameters:
tau = 1.99;  % recommended that tau is between 1 and 2
eta_mu = 1/4;   % how much to decrease mu at every step

%VERBOSE = false;  % no output
VERBOSE = 1;    % a little bit of output
VERBOSE = 2;    % even more output

if nargin < 6 || isempty(tol)
    tol = 1e-4;
end
if nargin < 5 || isempty(maxiter)
    maxiter = 500;
end
    
if length(n) == 1,
    n1 = n(1); n2 = n1;
elseif length(n) == 2,
    n1 = n(1); n2 = n(2);
end
%if n1*n2 < 100*100
    SMALLSCALE = true; 
    X = zeros(n1,n2);
%else
%    SMALLSCALE = false;
%end


m = length(Omega); [temp,indx] = sort(Omega); 
incre = 5;
r = 1; s = r + 1;  % estimate new rank
normb = norm(b);

[i, j] = ind2sub([n1,n2], Omega);
G = sparse(i,j,b,n1,n2,m);  % i.e. starting with X = 0;
mu = normest( G , 1e-2 );

% What the best way to multiply a sparse matrix?
[forwardType, transposeType] = findBestMultiply(G,.2);

U = zeros(n1,1);
V = zeros(n2,1);
S = 0;
relResid = 2;

if VERBOSE, fprintf('**************************************\n'); end
numiter = 0;
while mu > mu_final
    mu = max(mu * eta_mu,mu_final);
    if VERBOSE, fprintf('FPC, mu = %f\n',mu); end
    if VERBOSE == 1, fprintf('\tIteration:      '); end
    s = 2*r + 1;  % estimate new rank for next iteration
for k = 1:maxiter
    numiter = numiter + 1;
    if VERBOSE==1, fprintf('\b\b\b\b%4d',k);  end

    % Make routines for multiplying by a sparse matrix
    Gt = G';
    switch forwardType
        case 1, Gforward = @(x) G*x;
        case 2, Gforward = @(x) Gt'*x;
        case 3, Gforward = @(x) smvp(G,x);
    end
    switch transposeType 
        case 1, Gtranspose = @(x) Gt*x;
        case 2, Gtranspose = @(x) G'*x;
        case 3, Gtranspose = @(x) smvp(Gt,x);
    end
    
    % Y = X - tau*G
    Y = @(x) U*(S*(V'*x)) - tau*Gforward(x);
    Yt= @(x) V*(S*(U'*x)) - tau*Gtranspose(x);
    
    % Perform a SVD
    if SMALLSCALE
         [U,Sigma,V] = svd(full(X - tau*G),'econ');
         %[U,Sigma,V] = svdecon(full(X - tau*G));
    else
        OK = 0;
        while ~OK
            opts = []; 
            if ~isreal(G), opts.eta = 1e-16; end
%             opts.minSingValue = tau*mu;
%             opts.increaseK = 10;
%             [U,Sigma,V] = lansvd(Y,Yt,n1,n2,min(s+1,min(n1,n2)),'T',opts);
            [U,Sigma,V] = lansvd(Y,Yt,n1,n2,s,'L',opts);
            OK = (Sigma(s,s) <= tau*mu) || ( s == min(n1,n2) );
            s = min(2*s, min(n1,n2));
%             if ~OK, disp('increasing rank'); end
%             s = min(s + incre, min(n1,n2));
        end
    end
    % Shrink:
    sigma = diag(Sigma); r = sum(sigma > tau*mu); 
    U = U(:,1:r); V = V(:,1:r); 
    sigma = sigma(1:r) - tau*mu; 
    S = diag(sigma);
    s = r + 1;  % estimate new rank for next iteration
    
    % update P(X) and g(X) = P*(P(X)-b)
    if SMALLSCALE
        X = U*S*V'; x = X(Omega);
    else
        x = XonOmega(U*S,V,Omega);
    end
    resid = x - b;
    try
        updateSparse(G,resid,indx);
    catch
        l = lasterror;
        if strcmpi( l.identifier, 'MATLAB:UndefinedFunction')
            % mex file not installed, so do this instead:
            G = updateSparse_slow(G,resid,indx);
        else
            % some other error (unexpected)
            rethrow(lasterror)
        end
    end

    
    old_relResid = relResid;
    relResid = norm(resid)/normb;
    
    if VERBOSE == 2
        fprintf('iteration %4d, rank is %2d, rel. residual is %.1e\n',k,r,relResid);
    end
%     if (relResid < tol)
%         break
%     end
    
    % use stopping criteria 4.1 ("gtol")
    % need to compute spectral norm
    if ~rem(k,5)
        if SMALLSCALE
            sigma_max = norm( U*V' - full(G)/mu );
        else
            Y = @(x) U*(V'*x) - (G*x)/mu;
            Yt= @(x) V*(U'*x) - (G'*x)/mu;
            sigma_max = lansvd(Y,Yt,n1,n2,1,'L');
        end
        if sigma_max - 1 < tol
            break
        end
    end
    
    % Above stopping criteria rarely applies, so add this one also:
    if k > 1 && ( abs( relResid - old_relResid) / old_relResid ) < tol
        break
    end

    if (relResid > 1e5)
        disp('Divergence!');
        break
    end
    
end
if VERBOSE == 1
    fprintf('\n\tRelative Residual is %.3f%%\n',100*relResid); 
    fprintf('\tRank is %d\n',r);
end

end
