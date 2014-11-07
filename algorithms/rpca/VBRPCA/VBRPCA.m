function [X, A, B, E, varargout] = VBRPCA(Y, options)
% Variational Bayesian Robust Principal Component Analysis
% Author (a.k.a. person to blame): S. Derin Babacan
% Last updated: January 31, 2012
% Usage:
% [X, A, B, E, varargout] = VBRPCA(Y, options)
% Solves for A,B,E in
% Y = X + E + N
%   = AB' + E + N 
% where X, A, B are low rank and E is sparse, N is
% dense Gaussian noise
%
% If you're using this code, please acknowledge 
%    S. D. Babacan, M. Luessi, R. Molina, and A. K. Katsaggelos, 
%    "Sparse Bayesian Methods for Low-Rank Matrix Estimation," 
%    IEEE Transactions on Signal Processing, 2012.
%
% This code is not optimized for speed. 
% -----------------------------------------------------------------------
%                                INPUTS
% Y               : input matrix
% options         : All options are *optional*. The default values are set
%                   automatically.
%   
%   verbose       : output the progress? (0/1) default: 0. 
%   init          : Initialization method. 
%                   'rand': initialize A and B with random matrices
%                   'ml'  : Apply SVD to Y and initialize A and B using its factors.
%                   (default)
%   a_gamma0
%   b_gamma0      : hyperparameter values for gamma. default: 1e-6
%
%   a_alpha0
%   b_alpha0      : hyperparameter values for alpha. default: 0
%
%   a_beta0
%   b_beta0      : hyperparameter values for beta. default: 0
%   
%   inf_flag      : Flag for inference of components of E. default: 1
%                    1: Standard
%                    2: fixed point MacKay
%   MAXITER       : max number of iterations. default: 100
%   UPDATE_BETA   : update noise variance? (0/1). default: 1
%   DIM_RED       : prune irrelevant dimensions during iterations? (0/1).  default: 1
%   DIMRED_THR    : threshold to remove columns from A and B. Used only if
%                   DIM_RED is 1.  default: 1e4
%   initial_rank  : The initial rank of A and B to start. Default is full
%                   rank, calculated from the SVD of Y. A smaller value is
%                   helpful if the data dimensions is large.
%   X_true,E_true : Original X and E matrices for simulations. 
%                   If supplied and verbose=1,
%                   the errors are reported at each iteration
% -----------------------------------------------------------------------
%                                OUTPUTS
% X               : low rank component of Y
% A, B            : low rank factors in X = AB'
% E               : sparse component of Y
% optional outputs (see the paper for the full definitions):
%   varargout{1} = gamma    --- hyperparameters for A and B
%   varargout{2} = alpha    --- hyperparameters for E
%   varargout{3} = beta     --- noise inverse variance
%   varargout{4} = it       --- number of iterations
%   varargout{5} = t        --- running time
%   varargout{6} = Sigma_A  --- covariance matrix of A
%   varargout{7} = Sigma_B  --- covariance matrix of B
% -----------------------------------------------------------------------
%
% 
% Copyright (c): S. Derin Babacan, 2011
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.
%



%% Check input variables 
if ~exist('options','var')
    options = populate_vars([]);
else
    options = populate_vars(options);
end

verbose     = options.verbose;
a_gamma0    = options.a_gamma0;
b_gamma0    = options.b_gamma0;
a_beta0     = options.a_beta0;
b_beta0     = options.b_beta0;
a_alpha0    = options.a_alpha0;
b_alpha0    = options.b_alpha0;
MAXITER     = options.MAXITER;
inf_flag    = options.inf_flag;
UPDATE_BETA = options.UPDATE_BETA;
DIMRED      = options.DIMRED;

if DIMRED,
    if isfield(options, 'DIMRED_THR')
        DIMRED_THR = options.DIMRED_THR;
    else
        DIMRED_THR = 1e4;
    end
end

% For synthetic simulations
if isfield(options, 'X_true') & isfield(options, 'E_true')
    X_true = options.X_true;
    E_true = options.E_true;
end


%% Initialization
[m n] = size(Y);
mn = m*n;

Y2sum = sum(Y(:).^2);
scale2 = Y2sum / (mn);
scale = sqrt(scale2);

% Initialize A, B and E
switch options.init,
    % Maximum likelihood 
    case 'ml'
        
        [U, S, V] = svd(Y, 'econ');
        
        if strcmp(options.initial_rank, 'auto')
            r = min([m,n]);
        else
            r = options.initial_rank;
        end
        
        A = U(:,1:r)*(S(1:r,1:r)).^(0.5);
        B = (S(1:r,1:r)).^(0.5)*V(:,1:r)'; B = B';
        
        Sigma_A = scale*eye(r,r);
        Sigma_B = scale*eye(r,r);
        gammas = scale*ones(r,1);
        
        beta = 1./scale2;
        
        if UPDATE_BETA == 0 & isfield(options, 'beta') 
            beta = options.beta;
        end
        
        Sigma_E = scale*ones(r,r);
        alphas = ones(m,n)*scale;
        
    % Random initialization    
    case 'rand'
        r = min([m,n]);
        A = randn(m,r) * sqrt(scale);
        B = randn(n,r) * sqrt(scale);
        gammas = scale*ones(r,1);
        Sigma_A = eye(r) * scale;
        Sigma_B = eye(r) * scale;
        E = randn(m, n) * sqrt(scale);
        beta = 1./scale2;
        Sigma_E = scale*ones(r,r);
        alphas = ones(m,n)*scale;
        
end

X = A*B';
E = Y-X;

%% Iterations
tic
for it = 1:MAXITER,
    old_X = X;
    old_E = E;
    
    %% update X
    W = diag(gammas);
    %% A step
    
    if strcmp(options.mode, 'VB')
        Sigma_A = (beta*B'*B + W + beta*n*Sigma_B)^(-1);
        A = beta*(Y-E)*B*Sigma_A;
    elseif strcmp(options.mode, 'VB_app')
        Sigma_A = (beta*B'*B + W + beta*n*Sigma_B);
        A = (beta*(Y-E)*B ) / Sigma_A;
        Sigma_A = diag(1./diag(Sigma_A));
    elseif strcmp(options.mode, 'MAP')
        Sigma_A = (beta*B'*B + W + beta*n*Sigma_B);
        A = (beta*(Y-E)*B ) / Sigma_A;
        Sigma_A = zeros(size( Sigma_A ) ); 
    end
    
    %% B step
    if strcmp(options.mode, 'VB')
        Sigma_B = (beta*A'*A + W + beta*m*Sigma_A)^(-1);
        B = beta*(Y-E)'*A*Sigma_B;
    elseif strcmp(options.mode, 'VB_app')
        Sigma_B = (beta*A'*A + W + beta*m*Sigma_A);
        B = (beta*(Y-E)'*A)/Sigma_B;
        Sigma_B = diag(1./diag(Sigma_B));
    elseif strcmp(options.mode, 'MAP')
        Sigma_B = (beta*A'*A + W + beta*m*Sigma_A);
        B = (beta*(Y-E)'*A)/Sigma_B;
        Sigma_B = zeros(size( Sigma_B ) );
    end
    
    X = A*B';
   
    %% E-step
    Sigma_E = 1./(alphas+beta);
    E = beta*(Y-X).*Sigma_E;
    
    
    %% Estimate alphas
    %----------------------------------------------------------------------------------------
    % Choose inference method (fixed point Mackay or standard)
    if inf_flag == 1, 
        % Standard
        alphas = 1./(E.^2 + Sigma_E);
    elseif inf_flag == 2,
        % Mackay
        alphas = (1-alphas.*Sigma_E + a_alpha0)./ (E.^2 + eps + b_alpha0);
    end
    
    
    %----------------------------------------------------------------------------------------
    %% estimate gammas
	gammas = (m + n + a_gamma0)./( diag(B'*B) + diag(A'*A)+ m*diag(Sigma_A) + n*diag(Sigma_B)+ b_gamma0);

    % MacKay update for gamma? (not tested extensively)
%     gammas = ( m+n - m*diag(Sigma_A) - n*diag(Sigma_B) + a_gamma0 )./ (diag(B'*B) + diag(A'*A)+ b_gamma0);
    
    %% estimate beta? 
    if  UPDATE_BETA,
        err = sum(sum( abs(Y - X - E).^2 ) ) + n*trace(A'*A*Sigma_B) + m*trace(B'*B*Sigma_A) + m*n*trace(Sigma_A*Sigma_B) + sum(sum((Sigma_E)));
        beta = (m*n + a_beta0)/(err+b_beta0);
    end
    
    %% Prune irrelevant dimensions?
    if DIMRED,
        MAX_GAMMA = min(gammas) * DIMRED_THR;
        
        if sum(find(gammas > MAX_GAMMA)),
            indices = find(gammas <= MAX_GAMMA);
            
            A = A(:,indices);
            B = B(:,indices);
            gammas = gammas(indices);
            
            Sigma_A = Sigma_A(indices,indices);
            Sigma_B = Sigma_B(indices,indices);
            
            
            [m r] = size(A);
            [r n] = size(B'); % r is the current rank estimate
        end
        
    end
    
    
    %% Check convergence and display progress
    Xconv = norm(old_X-X,'fro')/norm(old_X,'fro');
    
    if verbose,
        
        Econv = norm(old_E - E,'fro')/ norm(old_E,'fro');
        
        if exist('X_true','var') % For synthetic simulations
            Xerr = norm( X - X_true , 'fro') / norm(  X_true , 'fro');
            Eerr = norm( E - E_true , 'fro') / norm(  E_true , 'fro');
            fprintf('it %d: X Error = %g, Xconv = %g, E Error = %g, Econv = %g, beta = %g, r = %d\n', it, Xerr, Xconv, Eerr, Econv, beta, r);
            
        else
            fprintf('it %d: Xconv = %g, Econv = %g, beta = %g, r = %d\n', it, Xconv, Econv, beta, r);
            
        end
    end
        
    % Check for convergence (do at least 5 iterations)
    if it > 5 & Xconv < options.thr 
        break;
    end
    
end
t = toc;

varargout{1} = gammas;
varargout{2} = alphas;
varargout{3} = beta;
varargout{4} = it;
varargout{5} = t;
varargout{6} = Sigma_A;
varargout{7} = Sigma_B;



%% Function to populate options.
function [options] = populate_vars(options)

if isempty(options),

    options.init         = 'ml';
    options.verbose      = 1;
    options.a_gamma0     = 1e-6;
    options.b_gamma0     = 1e-6;
    options.a_alpha0     = 0;
    options.b_alpha0     = 0;
    options.a_beta0      = 0;
    options.b_beta0      = 0;
    options.MAXITER      = 100;
    options.DIMRED       = 1;
    options.inf_flag     = 2;
    options.UPDATE_BETA  = 1;
    options.mode         = 'VB';
    options.thr          = 1e-5;
    options.initial_rank = 'auto';
    
else
    
    if ~isfield(options, 'init')
        options.init = 'ml';
    end
    
    if ~strcmp(options.init,'ml') | ~strcmp(options.init,'rand')
        options.init = 'ml';
    end
    
    if ~isfield(options, 'verbose')
        options.verbose   = 0;
    end
    
    if ~isfield(options, 'a_gamma0')
        options.a_gamma0   = 1e-6;
    end
    
    if ~isfield(options, 'b_gamma0')
        options.b_gamma0   = 1e-6;
    end
    
    if ~isfield(options, 'a_alpha0')
        options.a_alpha0   = 0;
    end
    
    if ~isfield(options, 'b_alpha0')
        options.b_alpha0   = 0;
    end
    
    if ~isfield(options, 'a_beta0')
        options.a_beta0   = 0;
    end
    
    if ~isfield(options, 'b_beta0')
        options.b_beta0   = 0;
    end
    
    if ~isfield(options, 'MAXITER')
        options.MAXITER   = 100;
    end
    
    if ~isfield(options, 'DIMRED')
        options.DIMRED   = 1;
    end
    
    if ~isfield(options, 'inf_flag')
        options.inf_flag    = 2;
    end
    
    if options.inf_flag ~= 1 | options.inf_flag ~= 2,
        options.inf_flag = 2;
    end
    
    if ~isfield(options, 'UPDATE_BETA')
        options.UPDATE_BETA = 1;
    end
    
    if ~isfield(options, 'mode')
        options.mode = 'VB';
    end
    
    if ~isfield(options, 'thr')
        options.thr = 1e-5;
    end
    
    if ~isfield(options, 'initial_rank')
        options.initial_rank = 'auto';
    end
    
end


