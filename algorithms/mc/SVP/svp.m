function [U S V t err] = svp (Omega, b, m, n, k, params)
%Return the result of solving min \| b - X(Omega) \|_F^2 , s.t rank(X) <= k
%Inputs:
%   Omega - indices of known entries.
%   b - values of known entries
%   k - approximate rank. We only work with matrices of rank up to rank k.
%   params - parameter structure containing the following fields:
%       tol,vtol,mxitr - stopping criteria
%           defaults = 1e-2 (tol), 1e-3 (vtol), 100 (mxitr)
%       verbosity - verbosity level
%           default = 1
%       eta - step size
%           default = 3/(4p) where p=fraction of entries available
%Outputs:
%   U - Left singular vectors of the optimal matrix X
%   S - Singular values of X
%   V - Right singular vectors of X
%   t - Number of iterations required
%   err - RMSE over the known values
%
%Reference:
% Raghu Meka, Prateek Jain, Inderjit S. Dhillon. 
% "Guaranteed Rank Minimization via Singular Value Projection"
% Arxiv:0909.5457
%Written by: Prateek Jain (pjain@cs.utexas.edu) and Raghu Meka (raghu@cs.utexas.edu)
%Last updated: October, 2009


if(~exist('Omega','var')||~exist('b','var')||~exist('m','var')||~exist('n','var')||~exist('k','var'))
    error('ErrorSVP:InvalidInputs','Not enough inputs.\n Usage: [U S V t err] = svp (Omega, b, m, n, k, params)\n Type "help svp" for more information');
end
if(~isreal(m)||~isreal(n)||~isreal(k))
    error('ErrorSVP:InvalidInputs',...
        'Incorrect input(s). m, n, k should all be real\nUsage: [U S V t err] = svp (Omega, b, m, n, k, params)\n Type "help svp" for more information');
end
if(length(Omega)~=length(b))
    error('ErrorSVP:InvalidInputs',...
        'Incorrect input. Omega and b have equal number of elements. \nUsage: [U S V t err] = svp (Omega, b, m, n, k, params)\n Type "help svp" for more information');
end
if(~isempty(find(Omega>m*n, 1)))
    error('ErrorSVP:InvalidInputs',...
        'Incorrect input. Some elements of Omega out of range (i.e., greater than m*n). \nUsage: [U S V t err] = svp (Omega, b, m, n, k, params)\n Type "help svp" for more information');
end

p=length(Omega)/(m*n); %sampling density

if(~exist('params','var'))
    params=struct();
end
params = SetDefaultParams(params,p);
tol=params.tol;
vtol=params.vtol;
maxtol=params.maxtol;
mxitr=params.mxitr;
verbosity=params.verbosity;
eta=params.eta; %Step Size



% Initialize X=U*S*V'=0
U = zeros(m,k);
S = zeros(k,1);
V = zeros(n,k);

% Compute indices for the known entries
[I, J] = ind2sub([m,n], Omega);

oerr=+Inf;
t = 1;

if verbosity == 1
    fprintf('\nIteration:   ');
end

while  t <= mxitr

    % Compute value of current iterate X=U*S*V' over the known entries
    X_Omega=compute_X_Omega(U*diag(S), V, Omega);

    % Compute RMSE and break if RMSE is small enough or not enough progress
    err=norm(X_Omega - b,'fro')/sqrt(length(b));

    if verbosity == 1
        fprintf('\b\b\b\b%4d',t);
    elseif verbosity >1
        fprintf('Iteration %d: Error = %f\n', t, err);
    end

    if err < tol || abs(oerr - err) < vtol
        break
    end

    if err-oerr >maxtol
      error(['Divergence!!! Current step size: ',num2str(eta),'. Decrease the step size eta (params.eta) and try' ...
	     ' again']);
    end
    oerr = err;

    % Compute a step in direction of gradient, i.e., Y=eta P_omega(M-X)
    y=eta*(b - X_Omega);

    % Take a step in the direction of gradient
    % and compute projection over top k singular vectors
    % i.e, U*S*V'=P_k(U*S*V'+Y)
    [U S V] = fastsvd(U, S, V, y, I, J, k);

    t = t+1;

end
fprintf('\n');

function [U,S,V] = fastsvd(U,S,V,y,I,J,k)
%Find the singular values of U S V' + Y
%Y is a sparse matrix
m = size(U,1);
n = size(V,1);
Afunc = @(x) (U*(S.*(V'*x)) + compOmegaYx(m,n,y,I,J,x));
Atfunc = @(x) (V*(S.*(U'*x)) + compOmegaYtx(m,n,y,I,J,x));
[U,Sigma,V] = lansvd(Afunc,Atfunc, m,n,k,'L');
S = diag(Sigma);