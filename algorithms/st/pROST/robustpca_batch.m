% Copyright (c) Technische Universität München, 2012
% All rights reserved.
%
% Author: Clemens Hage
% Contact: hage@tum.de
% Date: 2012-11-12
%
% Code based on the publication
% 'Robust PCA and subspace tracking from incomplete observations using L0-surrogates'
% by Clemens Hage and Martin Kleinsteuber
% Research Group for Geometric Optimization and Machine Learning
% Technische Universität München
% Preprint available at http://arxiv.org/abs/1210.0805
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


function [L,U,Y] = robustpca_batch(X,k,pentype,varargin)
%
% This function batch-processes a given data set and estimates the
% underlying subspace of dimension k. Can choose from three different 
% L0-surrogate cost functions.
% The algorithm approximately solves the problem min ||X - L||_0 s.t. 
% rank(L) <= k via minimization of h(X-UY), where U is an element of 
% St_(m,k), Y has dimensions k x n and h is a sparsifying penalty function.
% I alternating minimzation steps are performed. Firstly, the subspace 
% estimate is refined via Conjugate Gradient with backtracking line search
% over P=UU' on the Grassmannian (cost function f1), then Y is optimized 
% using CG in Euclidean space (cost function f2). The penalty function's
% smoothing parameter mu is reduced after each alternating minmization step.
%
% Mandatory inputs:
%
%       X - the input data of dimension m x n
%       k - the desired dimension of the underlying subspace. Must be an
%       integer value larger than 1.
%
%       pentype - the desired penalty function type. Select between the
%       following options:
%       'lpnorm' - a smoothed version of the Lp-norm:
%       X |-> Sum(j=1..n) Sum(i=1..m) (x_ij^2 + mu)^(p/2), 0 < p < 1
%       'lognorm' - a logarithmic penalty function:
%       X |-> Sum(j=1..n) Sum(i=1..m) log(1 + (x_ij^2 / mu))
%       'atansquare' - squared atan penalty function
%       X |-> Sum(j=1..n) Sum(i=1..m) atan^2(x_ij^2 / mu)
%       
%
%
% Optional parameters (default value):
%       I - the number of alternating iteration steps. (10)
%
%       mu_0 - initial value for the smoothing parameter. (depends on the
%       particular penalty function)
%
%       mu_I - final value for the smoothing parameter. (depends on the
%       particular penalty function)
%
%       p - p value for the Lp-norm (0.5)
%
%       Omega - The support set of the observation (all data observed)
%
%       L_0 - the ground truth of the underlying subspace (unknown)

%% Parse inputs
P = inputParser;

P.addRequired('X', @isnumeric);
P.addRequired('k', @(x)isnumeric(x) && x > 1 && x <= min(size(X)));
validpentypes={'lpnorm','lognorm','atansquare'};
P.addRequired('pentype',@(x)any(strcmp(x,validpentypes)));

switch(pentype)
    case('lpnorm')
        mu_0_default = 0.9;
        mu_I_default = 1e-8;
    case('lognorm')
        mu_0_default = 1;
        mu_I_default = 0.01;
    case('atansquare')
        mu_0_default = 1;
        mu_I_default = 0.1;
end

P.addOptional('I', 10, @isnumeric);
P.addOptional('mu_0',mu_0_default, @(x)isnumeric(x) && x > 0);
P.addOptional('mu_I',mu_I_default, @(x)isnumeric(x) && x > 0);
P.addOptional('p',0.5, @(x)isnumeric(x) && 0 < x && x < 1);
P.addOptional('Omega', [],@(x)islogical(x));
P.addOptional('L_0', [],@(x)isnumeric(x));

% parse the inputs
P.parse(X,k,pentype,varargin{:});

% initialize default parameter set (hard-coded parameters)
params=setDefaultParams();

% add all adjustable parameters

% the data set
params.X = X;
% dimension of the data samples
params.m = size(X,1);
% desired subspace dimension
params.k = k;
% number of alternating minimization runs
params.I = P.Results.I;
% type of l0-surrogate penalty function
params.penfun.type = pentype;
% initial and final value for the smoothing parameter
params.mu_0 = P.Results.mu_0;
params.mu_I = P.Results.mu_I;
% initialize smoothing factor
params.penfun.mu = params.mu_0;

% initialize penalty function parameters depending on the penalty function
switch(pentype)
    case('lpnorm')
        params.penfun.p = P.Results.p;
        params.penfun.mu = params.mu_0;
    case('lognorm')
        params.penfun.mu = params.mu_0;
    case('atansquare')
        params.penfun.mu = params.mu_0;
end

if(isempty(P.Results.L_0))
    params.GT_KNOWN = 0;
else
    params.GT_KNOWN = 1;
    L_0 = P.Results.L_0;
end

% check if obeservation is full or partial
params.FULL_OBS = isempty(P.Results.Omega);

if(~params.FULL_OBS)
    params.Omega = P.Results.Omega;
end

% shrinkage factor for mu
params.c_mu = (params.mu_I/params.mu_0)^(1/(params.I-1));

%% Modify hard-coded parameters. Refer to setDefaultParams() for parameter list
% params.VERBOSE = 1; % for more output
% params.PCAINIT = 0; % random initialization
% params.NONMONOTONE = 1; % nonmonotone line search


%% INITIALIZATION
% check if data is fully or partially known and initialize randomly or with
% PCA estimate.
if(params.FULL_OBS)
    if(params.PCAINIT)
        [U,~,~]=svds(X,k);
    else
        U=orth(randn(params.m,k));
    end
    Y = U' * X;
else
    X_Omega = full(params.Omega.*X);
    if(params.PCAINIT)
        [U,~,~]=svds(X_Omega,k);
    else
        U=orth(randn(params.m,k));
    end
    Y = U' * X_Omega;
end

% initial subspace estimate
L = U * Y;

%% ALTERNATING MINIMIZATION

for i = 1:params.I
    
    display(['Minimization pass #' num2str(i) ', mu =' num2str(params.penfun.mu)]);
    
    % optimize over P (via U)
    [U,it] = min_f1(U,L,params);
    
    if(params.VERBOSE)
        display(['Optimized P in ' num2str(it) ' iterations']);
    end
    
    % optimize over Y
    [Y,it] = min_f2(U,L,params);
    
    if(params.VERBOSE)
        display(['Optimized Y in ' num2str(it) ' iterations']);
    end
    
    % update L
    L = U*Y;
    
    if(params.GT_KNOWN)
        if(params.FULL_OBS)
            display(['relative error: ' num2str(norm(L-L_0,'fro') / norm(L_0,'fro'))])
        else
            display(['relative error: ' num2str(norm(params.Omega.*(L-L_0),'fro') / norm(params.Omega.*L_0,'fro'))])
        end
    end
    
    % shrink smoothing parameter
    params.penfun.mu = params.c_mu * params.penfun.mu;
end

end

%% MINIMIZING U (f1)
function fval = f1(U,L,params)
% evaluates cost function f1

if(params.FULL_OBS)
    fval = penfun(params.X - U*(U'*L),params.penfun);
else
    fval = penfun(params.Omega.*(params.X - U*(U'*L)),params.penfun);
end
end

function grad = delta_f1(U,L,params)
% computes the symmetrized gradient of f1

if(params.FULL_OBS)
    deltah = delta_penfun(params.X - U*(U'*L),params.penfun);
else
    deltah = delta_penfun(params.Omega.*(params.X - U*(U'*L)),params.penfun);
end
grad = - L*deltah';

% symmetrize
grad = 1/2 * (grad + grad');
end

function [U_new,it] = min_f1(U,L,params)
m = params.m;
k = params.k;

FORCE_INIT=0;

% initialize history to store cost function progress
f_history = zeros(params.f1.maxit,1);

V = qr_positive(U); % so that V' * P * V is the standard projector
U_orth = V(:,(k+1):m); % orthogonal complement of U

% G is the lower left part of the matrix obtained by conjugating the
% Riemannian gradient Gamma = pi_P(grad_f1) with V
G = U_orth' * (delta_f1(U,L,params) * U);

H = -G; %initial search direction

GH_old = innerprod(G,H);
t = params.f1.t_init; % step-size initialization

% moving average over function values
f_av = 0;

for it = 1:params.f1.maxit
    
    [Theta_H,R] = qr_positive(H);
    R = R(1:min(k,m-k),:); % R is upper triangular kxk-matrix
    
    % need only the first 2k cols of Theta to compute U
    Theta = V * [eye(k) zeros(k); zeros(m-k,k) Theta_H(:,1:k)];
    
    % compute initial stepsize
    if(FORCE_INIT)
        FORCE_INIT=0;
    else
        t_init = min(10 * t * abs(GH_old / innerprod(G,H)), params.f1.t_init);
        if(mod(it,params.f1.reset_rate)==0)
            t_init= params.f1.t_init;
        end
    end
    
    t = bt_f1(U,L,G,H,Theta,R,t_init,params,f_av); % determine stepsize by backtracking line-search
    
    M = Mmatrix(t*R,k);
    Theta_M = qr_positive(M);
    % need only the first k cols for the computation of U
    Theta_M = Theta_M(:,1:k);
    
    % backup U
    U_old = U;
    
    % update U
    U = Theta * Theta_M;
    
    % evaluate the cost function
    f_history(it) = f1(U,L,params);
    
    % CONVERGENCE CHECK
    if(it > params.f1.bufsize) % fill up the cost function buffer first
        % compute average over cost function buffer
        f_av = mean(f_history(it-params.f1.bufsize : it-1));
        
        % compute relative progress in the cost function
        diff = abs((f_av - f_history(it)) / f_av);
        
        if(params.VERBOSE)
            display(['relative progress in f1: ' num2str(diff)]);
        end
        
        % Ensure algorithm cannot increase cost function (if a too large
        % stepsize is chosen erroneously)
        if(diff > 1)
            U = U_old;
            t_init = t_init / 10;
            FORCE_INIT=1;
            continue;
        end
        
        if(diff <= params.f1.thresh || t <= params.f1.t_min)
            % algorithm converged if relative progress or the step size 
            % fall below certain threshold
            break;
        end
    end
    
    % needed for backtracking initialization
    GH_old = innerprod(G,H);
    
    % transport G and H
    tauG = Theta_H * G;
    tauH = Theta_H * H;
    
    % update V, U_orth and G
    V = qr_positive(U);
    U_orth = V(:,(k+1):m);
    G = U_orth' * (delta_f1(U,L,params) * U);
    
    % Hestenes-Stiefel formula
    gamma = innerprod((G-tauG),G) / (innerprod((G-tauG),tauH) + eps);
    
    % update search direction
    H = -G + gamma * tauH;
    
    % Reset CG search direction after (m-k)*k iterations
    if(mod(it+1,(m-k)*k)==0)
        H = -G;
    end
end

U_new = U;

end

function t = bt_f1(U,L,G,H,Theta,R,t_init,params,f_av)
k = params.k;

% evaluate cost function at current point (t=0)
f = f1(U,L,params);

% set initial step size (will be multiplied with rho before first use)
t = t_init/params.f1.rho;

while (t > params.f1.t_min)
    t = params.f1.rho * t;
    
    % compute Theta_M for current t
    Theta_M = qr_positive(Mmatrix(t*R,k));
    
    % compute U for current t (Theta does not change)
    U_t = Theta * Theta_M(:,1:k);
    
    % evaluate cost function for current t
    f_t = f1(U_t,L,params);
    
    % check for break condition (Armijo rule)
    if(params.NONMONOTONE && f_av)
        if(f_t <= f_av + params.f1.c*t*innerprod(G,H))
            break;
        end
    else
        if(f_t <= f + params.f1.c*t*innerprod(G,H))
            break;
        end
    end
    
end
end

%% MINIMIZING Y (f2)
function fval = f2(U,Y,params)
% evaluates the cost function f2

if(params.FULL_OBS)
    fval = penfun(params.X - U*Y,params.penfun);
else
    fval = penfun(params.Omega.*(params.X - U*Y),params.penfun);
end
end

function grad = delta_f2(U,Y,params)
% computes the gradient of f2

if(params.FULL_OBS)
    grad = -U' * delta_penfun(params.X - U*Y,params.penfun);
else
    grad = -U' * delta_penfun(params.Omega.*(params.X - U*Y),params.penfun);
end
end

function [Y_new,it] = min_f2(U,L,params)
m = params.m;
k = params.k;

% flag is set to re-init the step size if algorithm behaves irregularly
FORCE_INIT=0;

% initialize history to store cost function progress
f_history = zeros(params.f2.maxit,1);

% initialize Y = U_i+1' * U_i* Y_i
Y = U' * L;

G = delta_f2(U,Y,params);
H = -G;

GH_old = innerprod(G,H);
t = params.f2.t_init; % step-size initialization

for it = 1:params.f2.maxit
    
    % compute initial stepsize
    if(FORCE_INIT)
        FORCE_INIT=0;
    else
        % ensure initial step size does not exceed max value
        t_init = min(10 * t * abs(GH_old / innerprod(G,H)), params.f2.t_init);
        if(mod(it,params.f2.reset_rate)==0)
            t_init= params.f2.t_init;
        end
    end
    
    t = bt_f2(U,Y,G,H,t_init,params); % determine stepsize by backtracking line-search
    
    % backup Y
    Y_old = Y;
    
    % update Y
    Y = Y + t * H;
    
    % evaluate the cost function
    f_history(it) = f2(U,Y,params);
    
    % CONVERGENCE CHECK
    if(it > params.f2.bufsize) % fill up the cost function buffer first
        % compute average over cost function buffer
        f_av = mean(f_history(it-params.f2.bufsize : it-1));
        
        % compute relative progress in the cost function
        diff = abs((f_av - f_history(it)) / f_av);
        
        if(params.VERBOSE)
            display(['relative progress in f2: ' num2str(diff)]);
        end
        
        % Ensure algorithm cannot increase cost function (may happen if a
        % large stepsize is chosen erroneously)
        if(diff >1)
            Y = Y_old;
            t_init = t_init / 10;
            FORCE_INIT=1;
            continue;
        end
        
        if(diff <= params.f2.thresh || t <= params.f2.t_min)
            % algorithm converged if step size or relative progress fall 
            % below certain threshold
            break;
        end
    end
    
    % needed for backtracking initialization
    GH_old = innerprod(G,H);
    
    % backup G_i before updating G
    G_old = G;
    
    G = delta_f2(U,Y,params);
    
    % Hestenes-Stiefel formula
    gamma = innerprod(G-G_old,G) / (innerprod(G-G_old,H) + eps);
    
    % update search direction
    H = -G + gamma * H;
    
    % Reset CG search direction after (m-k)*k iterations
    if(mod(it+1,(m-k)*k)==0)
        H = -G;
    end
end

Y_new = Y;

end

function t = bt_f2(U,Y,G,H,t_init,params)

% set initial step size (will be multiplied with rho before first use)
t = t_init / params.f2.rho;

% evaluate cost function at current point (t=0)
f = f2(U,Y,params);

while(t > params.f2.t_min)
    
    t = params.f2.rho * t;
    
    % compute Y for current t
    Y_t = Y + t * H;
    % evaluate cost function for current t
    f_t = f2(U,Y_t,params);
    
    % Armijo rule
    if(f_t <= f + params.f2.c*t*innerprod(G,H))
        break;
    end
end
end

%% L0-SURROGATES AND THEIR GRADIENTS

function fval = penfun(A,pfparams)

switch pfparams.type
    case('lpnorm')
        fval = sum(sum(exp(pfparams.p/2 * log(A.*A + pfparams.mu))));
    case('lognorm')
        mu_inv=1/pfparams.mu;
        fval = sum(sum(log(1 + mu_inv*A.*A)));
    case('atansquare')
        mu_inv=1/pfparams.mu;
        fval = sum(sum(atan(mu_inv*A).^2));
end

end

function grad = delta_penfun(A,pfparams)

switch pfparams.type
    case('lpnorm')
        grad = pfparams.p*A.*exp((pfparams.p/2-1) * log(A.*A + pfparams.mu));
    case('lognorm')
        mu_inv=1/pfparams.mu;
        grad = 2*mu_inv*A ./ (1+mu_inv*A.*A);
    case('atansquare')
        mu_inv=1/pfparams.mu;
        grad = 2*mu_inv*atan(mu_inv*A) ./ (1 + mu_inv^2*A.*A);
end

end

%% HELPER FUNCTIONS
function [Q,R] = qr_positive(A)
% QR decomposition enforcing positive entries on the main diagonal of R

[Q,R] = qr(A);

dvec=sign(diag(R));
d=numel(dvec);

R(1:d,:)=bsxfun(@times,R(1:d,:),dvec);
Q(:,1:d)=bsxfun(@times,Q(:,1:d),dvec');

end

function M = Mmatrix(A,k)

M = [eye(k) -A'; A eye(k)];

end

function fval = innerprod(A,B)
% computes inner product tr(B'*A) without performing the full matrix multiplication
fval = reshape(B,numel(B),1)' * reshape(A,numel(A),1);

end

function params=setDefaultParams()
% This function initializes the parameter struct. Some parameters are 
% modified by the input parser

% FLAGS
% verbose mode
params.VERBOSE = 0;
% if data is fully observed
params.FULL_OBS = 1;
% if the ground truth is available
params.GT_KNOWN = 0;
% if common PCA should serve as an initial estimate
params.PCAINIT = 1;
% use nonmonotone line search
params.NONMONOTONE=0;

% BACKTRACKING PARAMETERS
% Minimization of f1
% Step size shrinkage rate
params.f1.rho = 0.5;
% Armijo rule variable
params.f1.c = 1e-6;
% stopping threshold
params.f1.thresh = 1e-4;
% max #iterations
params.f1.maxit = 100;
% function value buffer size (also minimum #performed iterations)
params.f1.bufsize = 3;
% step size initialization
params.f1.t_init = 1e-3;
% t_init reset rate
params.f1.reset_rate = 10; % start with t_init every 10 iterations
% minimum allowed step size
params.f1.t_min = 1e-16;

% Minimization of f2
% Step size shrinkage rate
params.f2.rho = 0.5;
% Armijo rule variable
params.f2.c = 1e-4;
% stopping threshold
params.f2.thresh = 1e-4;
% max #iterations
params.f2.maxit = 100;
% function value buffer size (also minimum #performed iterations)
params.f2.bufsize = 3;
% step size initialization
params.f2.t_init = 1;
% t_init reset rate
params.f2.reset_rate = 10; % start with t_init every 10 iterations
% minimum allowed step size
params.f2.t_min = 1e-16;
end