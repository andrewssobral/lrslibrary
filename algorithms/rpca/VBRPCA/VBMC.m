function [X, A, B, t, varargout] = VBMC(P, Y, options)
% Variational Bayesian Low Rank Matrix Completion
% Author (a.k.a. person to blame): S. Derin Babacan
% 
% Last updated: January 31, 2012
% Usage:
% [X, A, B, t, varargout] = VSBL_MAT(P, Y, options)
% Solves for X = AB' in
% Y = P(X) + N = P(AB') + N  where X, A, B are low rank, 
% N is dense Gaussian noise, and P is a binary sampling matrix.
%
% If you're using this code, please acknowledge 
% Reference: 
%    S. D. Babacan, M. Luessi, R. Molina, and A. K. Katsaggelos, 
%    "Sparse Bayesian Methods for Low-Rank Matrix Estimation," 
%    IEEE Transactions on Signal Processing, 2012.
% 
% This code is not optimized for speed. It implements the full variational
% Bayesian inference described in the paper, without any manipulation of
% the matrices. These are described in the paper above but not implemented
% in this version of the code. These manipulations might lead to significant 
% decrease in running times depending on the size of matrix Y. 
% 
% -----------------------------------------------------------------------
%                                INPUTS
%   P               : binary (0-1) sampling matrix same size as Y
%   Y               : input matrix
%   options         : All options are *optional*. The default values are set
%                     automatically.
%   
%   verbose       : output the progress? (0/1) default: 0. 
%   init          : Initialization method. 
%                   'rand': initialize A and B with random matrices
%                   'ml'  : Apply SVD to Y and initialize A and B using its factors.
%                   (default)
%   a_gamma0
%   b_gamma0      : hyperparameter values for gamma. default: 1e-6
%
%
%   a_beta0
%   b_beta0      : hyperparameter values for beta. default: 0
%   
%   inf_flag      : Flag for inference of components of E. default: 1
%                    1: Standard
%                    2: fixed point MacKay
%   MAXITER       : max number of iterations. default: 100
%   UPDATE_BETA   : update noise variance? (0/1). default: 1
%   UPDATE_BETA_START : iteration number to start updating the noise variance
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
% X               : low rank component
% A, B            : low rank factors of X = AB'
% optional outputs (see the paper for the full definitions):
%   varargout{1} = gamma    --- hyperparameters for A and B
%   varargout{2} = beta     --- noise inverse variance
%   varargout{3} = it       --- number of iterations
%   varargout{4} = t        --- running time
%   varargout{5} = Sigma_A  --- covariance matrix of A
%   varargout{6} = Sigma_B  --- covariance matrix of B
% -----------------------------------------------------------------------
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

%% Check input variables, 
if ~exist('options','var')
    options = populate_vars([]);
else
    options = populate_vars(options);
end

verbose      = options.verbose;
a_gamma0     = options.a_gamma0;
b_gamma0     = options.b_gamma0;
a_beta0      = options.a_beta0;
b_beta0      = options.b_beta0;
MAXITER      = options.MAXITER;
inf_flag     = options.inf_flag;
UPDATE_BETA  = options.UPDATE_BETA;
DIMRED       = options.DIMRED;

if DIMRED,
    if isfield(options, 'DIMRED_THR')
        DIMRED_THR = options.DIMRED_THR;
    else
        DIMRED_THR = 1e4;
    end
end

UPDATE_BETA_START = options.UPDATE_BETA_START;

  
% For synthetic simulations
if isfield(options, 'X_true') 
    X_true = options.X_true;
end



%% Initialization
[m n] = size(Y);
mn = m*n;
pmn = length(find(P == 1));
p = pmn/mn;

Y2sum = sum(Y(:).^2);
scale2 = Y2sum / (mn);
scale = sqrt(scale2);

% Initialize A and B 
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
        
        Sigma_A = repmat( scale*eye(r,r), [1 1 m] );
        Sigma_B = repmat( scale*eye(r,r), [1 1 n] );

        gammas = 1./diag(S.^2);
        gammas = gammas(1:r);
        gammas = ones(r,1)*scale;
        
        ssscale2 = 1;
        beta = ssscale2./scale2;
        
        if UPDATE_BETA == 0 & isfield(options, 'beta') 
            beta = options.beta;
        end
        
    % Random initialization    
    case 'rand'
        
        if strcmp(options.initial_rank, 'auto')
            r = min([m,n]);
        else
            r = options.initial_rank;
        end
        
        A = randn(m,r) * sqrt(scale);
        B = randn(n,r) * sqrt(scale);
        gammas = scale*ones(r,1);
        
        Sigma_A = repmat( scale*eye(r,r), [1 1 m] );
        Sigma_B = repmat( scale*eye(r,r), [1 1 n] );
        
        if UPDATE_BETA == 0 & isfield(options, 'beta') 
            beta = options.beta;
        else
            beta = 1./scale2;
        end
        
                
end

X = A*B';


%% Iterations
tic
for it=1:MAXITER,
        
    old_X = X;

    Aw = diag(gammas);
    
    %% A step
    diagsA = 0;
    for i=1:m, %iterate over rows
        
        observed = find(P(i,:));
        Bi = B(observed,:);
        Sigma_A(:,:,i) = (beta*Bi'*Bi + beta*sum(Sigma_B(:,:,observed),3) + Aw)^(-1);
        A(i,:) = beta*Y(i,observed)*Bi*Sigma_A(:,:,i);
        diagsA = diagsA + diag( Sigma_A(:,:,i) );
        
    end
    
    
    %% B step
    diagsB = 0;
    for j=1:n, %Iterate over cols
        
        observed = find(P(:,j));
        Aj = A(observed,:);
        Sigma_B(:,:,j) = (beta*Aj'*Aj + beta*sum(Sigma_A(:,:,observed),3) + Aw)^(-1);
        B(j,:) = beta*Y(observed,j)'*Aj*Sigma_B(:,:,j);
        diagsB = diagsB + diag( Sigma_B(:,:,j) );
        
    end
    
    %% update X
    X = A*B';
    
    %% estimate gammas
    % Choose inference method (fixed point Mackay or standard)
    if inf_flag == 1, %Standard
        gammas = (m + n + a_gamma0)./( diag(B'*B) + diag(sum(Sigma_B,3)) + diag(A'*A)+ diag(sum(Sigma_A,3))+ b_gamma0);
    elseif inf_flag == 2, %Mackay
        gammas = (m + n + a_gamma0 - gammas.*( diag(sum(Sigma_A,3)) + diag(sum(Sigma_B,3)) ))./( diag(A'*A) + diag(B'*B)+b_gamma0);
    end
    
    
    %% estimate beta
    if  UPDATE_BETA & it>=UPDATE_BETA_START,
        

        err = sum(sum( abs(Y - P.*X).^2 ) );
        
        for l = 1: m
            observed = find(P(l, :));
            err = err + trace(B(observed,:)'*B(observed,:)*Sigma_A(:, :, l)) ...
                + trace(A(l, :)'*A(l, :) *sum(Sigma_B(:, :, observed), 3)) ...
                + trace( Sigma_A(:, :, l)*sum(Sigma_B(:, :, observed), 3));
            
        end
        
        beta = (pmn + a_beta0)/(err+b_beta0);
    end
    
    %% Prune irrelevant dimensions?   
    if DIMRED,
        MAX_GAMMA = min(gammas) * DIMRED_THR;
        
        if sum(find(gammas > MAX_GAMMA)),
            
            indices = find(gammas <= MAX_GAMMA);
            
            A = A(:,indices);
            B = B(:,indices);
            gammas = gammas(indices);
            
            Sigma_A = Sigma_A(indices,indices,:);
            Sigma_B = Sigma_B(indices,indices,:);
            
            [m r] = size(A);
            [n r] = size(B);
        end
        
    end
    
    %% Check convergence and display progress
    conv = sum(sum( ( old_X - X ).^2 ) ) / mn;

    if verbose, 
        
        if exist('X_true','var')
            % For synthetic simulations
            err = sum(sum( ( X - X_true ).^2 ) ) / sum(sum( ( X_true ).^2 ) );
            rms = sqrt( sum(sum( ( X - X_true ).^2 ) )/mn );
            relrec = norm(X_true-X,'fro')/norm(X_true,'fro');
            
            fprintf('it %d: Error = %g, beta = %g, conv = %g, r = %d\n', it, relrec, beta, conv, r);
        
        else
            
            fprintf('it %d: beta = %g, conv = %g, r = %d\n', it, beta, conv, r);
            
        end
    end
    
    % Check for convergence (do at least 5 iterations)
    if it > 5 & conv < options.thr 
        break;
    end
    
    
end
t = toc;

varargout{1} = gammas;
varargout{2} = beta;
varargout{3} = it;
varargout{4} = t;
varargout{5} = Sigma_A;
varargout{6} = Sigma_B;


%% Function to populate options.
function [options] = populate_vars(options)


if isempty(options)
    options.init         = 'ml';
    options.verbose      = 1;
    options.a_gamma0     = 1e-6;
    options.b_gamma0     = 1e-6;
    options.a_beta0      = 0;
    options.b_beta0      = 0;
    options.MAXITER      = 100;
    options.DIMRED       = 1;
    options.UPDATE_BETA  = 1;
    options.mode         = 'VB';
    options.thr          = 1e-7;
    options.initial_rank = 'auto';
    options.inf_flag     = 1;
    options.UPDATE_BETA_START = 1;
else
    
    if ~isfield(options, 'init')
        options.init = 'ml';
    end
    
    if ~strcmp(options.init,'ml') | ~strcmp(options.init,'rand')
        options.init = 'ml';
    end
    
    if ~isfield(options, 'verbose')
        options.verbose   = 1;
    end
    
    if ~isfield(options, 'a_gamma0')
        options.a_gamma0   = 1e-6;
    end
    
    if ~isfield(options, 'b_gamma0')
        options.b_gamma0   = 1e-6;
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
        options.inf_flag   = 1;
    end
    
    if options.inf_flag ~= 1 & options.inf_flag ~= 2,
        options.inf_flag = 1; % Standard inference
    end
    
    if ~isfield(options, 'UPDATE_BETA')
        options.UPDATE_BETA = 1;
    end
    
    if ~isfield(options, 'UPDATE_BETA_START') 
        options.UPDATE_BETA_START = 1;
    end

    
    if ~isfield(options, 'mode')
        options.mode = 'VB';
    end
    
    if ~isfield(options, 'thr')
        options.thr = 1e-7;
    end
    
    if ~isfield(options, 'initial_rank')
        options.initial_rank = 'auto';
    end
    
end



