%%  Copyright (C) 2011-2012 by Jun He, Laura Balzano, and Arthur Szlam
%  
%  This file is part of the GRASTA library.
%  It is provided without any warranty of fitness
%  for any purpose. You can redistribute this file
%  and/or modify it under the terms of the GNU
%  Lesser General Public License (LGPL) as published
%  by the Free Software Foundation, either version 3
%  of the License or (at your option) any later version.
%  (see http://www.opensource.org/licenses for more info)

function [ Unew, STATUSnew, OPTSnew ] = grasta_stream( y_Omega, idx, U0, STATUS,OPTIONS, OPTS)
%  grasta_stream is the streaming version of GRASTA which can be used for
%  online data processing
%
%   [1] Jun He, Laura Balzano, and John C.S. Lui. Online Robust Subspace Tracking from Partial Information
%       http://arxiv.org/abs/1109.3827
%   [2] Jun He, Laura Balzano, and Arthur Szlam. Incremental gradient on the grassmannian for online foreground 
%       and background separation in subsampled video. In IEEE Conference on Computer Vision and Pattern Recognition 
%       (CVPR), June 2012.
%
%  Usage:
%       [ Unew, STATUSnew, OPTSnew ] = grasta_stream( y_Omega, idx, U0, STATUS, OPTIONS, OPTS)
%   Inputs:
%       y_Omega: current partial observation of a full data vector
%       idx:     the observed indices of y_Omega
%       U0:      the current estimated subspace
%       STATUS:  the current status of GRASTA
%           last_mu: \mu_{t-1} for adaptive step rule
%           step_scale: the estimated step_scale constant for adaptive
%       step-size rule
%           lelve: the current level of the "Multi-Level" adaptive rule
%           last_gamma and last_w: previous gradient = last_gamma*last_w'
%
%       OPTIONS: the options of GRASTA    
%           CONSTANT_STEP: default 0, will use multi-level step-size rule; > 0 will
%           use constant step-size rule
%           DIM_M : the ambient dimension
%           RANK  : the estimated low-rank of the underlying system
%           rho : ADMM constant step 
%           ITER_MAX: the max_iter for ADMM solver
%           TOL : the stopping tolerance of admm_srp / mex_srp
%
%       OPTS:    the current options for the subproblem of GRASTA
%           refer to mex_srp.cpp or admm_srp.m
%   Outputs:
%       Unew:    the updated subspace 
%       STATUSnew: the updated running status
%       OPTSnew: the updated options for the subproblem of GRASTA
%
%  Author: Jun He, Laura Balzano
%  Email:  hejun.zz@gmail.com, sunbeam@ece.wisc.edu
%  Date:   Sept. 01, 2012
%
%

LEVEL_FACTOR        = 2;

if isfield(OPTIONS,'MIN_MU'),
    MIN_MU = OPTIONS.MIN_MU;
else
    MIN_MU = 1;
end

if isfield(OPTIONS,'MAX_MU'),
    MAX_MU = OPTIONS.MAX_MU;
else
    MAX_MU = 15;
end


if isfield(OPTIONS,'USE_MEX'),
    USE_MEX = OPTIONS.USE_MEX;
else
    USE_MEX = 1;
end

if isfield(OPTIONS,'CONSTANT_STEP'),
    CONSTANT_STEP = OPTIONS.CONSTANT_STEP;
else
    CONSTANT_STEP = 0;
end


if isfield(OPTIONS,'DIM_M'),
    DIM_M = OPTIONS.DIM_M;
else
    error('Should specify OPTIONS.DIM_M data ambient dimension!!!\n');
end

if isfield(OPTIONS,'ITER_MIN'),
    MIN_ITER = OPTIONS.ITER_MIN;
else
    MIN_ITER = 5;
end

if isfield(OPTIONS,'ITER_MAX'),
    ITER_MAX = OPTIONS.ITER_MAX;
else
    ITER_MAX = 60;
end

if isfield(OPTIONS,'TOL'),
    TOL = OPTIONS.TOL;
else
    TOL = 1e-6;
end


if isfield(OPTIONS,'MAX_LEVEL'),
    MAX_LEVEL = OPTIONS.MAX_LEVEL;
else
    MAX_LEVEL = 20;
end



if isfield(OPTIONS,'RANK'),
    RANK = OPTIONS.RANK;
else
    error('Should specify OPTIONS.RANK!!!\n');
end

if isfield(OPTIONS,'QUIET'),
    QUIET = OPTIONS.QUIET;
else
    QUIET = 0;
end

DEFAULT_MU_HIGH     = (MAX_MU-1)/2; 
DEFAULT_MU_LOW      = MIN_MU + 2;


% If we first call GRASTA_stream then do the following initilization 
if STATUS.init == 0,    
    STATUS.init         = 1; % Do not enter this initial part any longer   
    STATUS.curr_iter    = 0; % For debug
    
    STATUS.last_mu      = MIN_MU;        
    STATUS.level        = 0;
    STATUS.step_scale   = 0;
    STATUS.last_w       = zeros(RANK,1);
    STATUS.last_gamma   = zeros(DIM_M,1);

    OPTS.TOL        = TOL;
    OPTS.MAX_ITER   = MIN_ITER;  % the max iteration of ADMM at level=0
    OPTS.QUIET      = 1;    

    if isfield(OPTIONS,'rho'),
        OPTS.RHO = OPTIONS.rho;
    else
        OPTS.RHO = 1.8;
    end
    
    U0 = orth(randn(DIM_M,RANK));
end

%%%%%%%%%%%%%%%
% main framework of GRASTA
U_Omega = U0(idx,:);

if USE_MEX,
    [s_t, w, ldual] = mex_srp(U_Omega, y_Omega, OPTS);
else
    [s_t, w, ldual] = admm_srp(U_Omega, y_Omega, OPTS);
end

gamma_1 = ldual;% + OPTS.RHO*(U_Omega*w + s_t - y_Omega);
UtDual_omega = U_Omega' * gamma_1;
gamma_2 = U0 * UtDual_omega;
gamma = zeros(DIM_M,1);
gamma(idx) = gamma_1;
gamma = gamma - gamma_2;

gamma_norm = norm(gamma);
w_norm     = norm(w);
sG = gamma_norm * w_norm;


% Here we use the adaptive step-size rule for SGD. The adaptive can work
% well for both dynamic and static subspace tracking tasks

% 1. determine the step scale from the first observation
if ~STATUS.step_scale,
    STATUS.step_scale = 0.5*pi*(1+MIN_MU)/sG;
    
    if ~QUIET,
        fprintf('Level 0: %.2e\n',STATUS.step_scale);
    end
end

% 2. inner product of previous grad and current grad
grad_ip = trace(STATUS.last_w * (STATUS.last_gamma' * gamma) * w');

%%% avoid inner product too large
normalization = norm(STATUS.last_gamma * STATUS.last_w','fro') * norm(gamma * w','fro');
if normalization == 0,
    grad_ip_normalization = 0;
else
    grad_ip_normalization = grad_ip/normalization;
end
%%%


% 3. if the two consecutive grad in the same direction, we take a larger
% step along the gradient direction, otherwise take a small step along the
% gradient direction
STATUS.last_mu = max(STATUS.last_mu + sigmoid(-grad_ip_normalization) , MIN_MU);

if CONSTANT_STEP > 0,
    t = CONSTANT_STEP;
else
    % should not take a step larger than pi/2
    t = STATUS.step_scale * LEVEL_FACTOR^(-STATUS.level) * sG / (1+STATUS.last_mu);
    if t>=pi/3,
        t = pi/3;
    end
    
    % Adjust the level
    bShrUpd = 0;
    if STATUS.last_mu <= MIN_MU,
        if STATUS.level > 1,
            bShrUpd = 1;
            STATUS.level = STATUS.level - 1;
            
            if ~QUIET,
                fprintf('multi-level adaption - decreasing, t:%.2e, vectors: %d, level: %d\n',...
                    t,STATUS.curr_iter, STATUS.level);
            end
            STATUS.curr_iter  = 0;
        end
        
        STATUS.last_mu    = DEFAULT_MU_LOW;
    elseif STATUS.last_mu > MAX_MU,
        if STATUS.level < MAX_LEVEL,
            bShrUpd = 1;
            STATUS.level = STATUS.level + 1;
            
            if ~QUIET,
                fprintf('multi-level adaption - increasing, t:%.2e, vectors: %d, level: %d\n',...
                    t, STATUS.curr_iter,  STATUS.level);
            end
            STATUS.curr_iter  = 0;
            STATUS.last_mu    = DEFAULT_MU_HIGH;
        else
            STATUS.last_mu    = MAX_MU;
        end
    end
    
    if bShrUpd,
        if STATUS.level>=0 && STATUS.level <4,          % [0,4)
            OPTS.MAX_ITER = MIN_ITER;
        elseif STATUS.level>=4 && STATUS.level <7,      % [4,7)
            OPTS.MAX_ITER = min(MIN_ITER*2, ITER_MAX);
        elseif STATUS.level>=7 && STATUS.level <10,     % [7,10)
            OPTS.MAX_ITER = min(MIN_ITER*4, ITER_MAX);
        elseif STATUS.level>=10 && STATUS.level <14,    % [10,14)
            OPTS.MAX_ITER = min(MIN_ITER*8, ITER_MAX);
        else
            OPTS.MAX_ITER = ITER_MAX;               % [14,...)
        end
        
        if ~QUIET,
            fprintf('Will use %d ADMM iterations in level %d\n',OPTS.MAX_ITER, STATUS.level);
        end
    end
    
end

% 4. update the gradient for further step size update
STATUS.last_gamma  = gamma;
STATUS.last_w      = w;

STATUS.grad_ip = grad_ip_normalization; % just for debugging


% Take the gradient step along Grassmannian geodesic.
alpha = w/w_norm;
beta  = gamma/gamma_norm;
step = (cos(t)-1)*U0*(alpha*alpha')  - sin(t)*beta*alpha';

U0 = U0 + step;

%%

STATUS.s_t   = s_t;
STATUS.w     = w;
STATUS.ldual = ldual;
STATUS.SCALE = 1;
STATUS.curr_iter = STATUS.curr_iter + 1;

STATUS.grasta_t = t;

%
%%%%%%%%%%%%%%%%%%%%%%

Unew = U0;
STATUSnew = STATUS;
OPTSnew = OPTS;
end


%% Function of Sigmoid

function fval = sigmoid(x)
FMIN = -1; FMAX = 1;
omega = 0.1;

fval = FMIN + (FMAX - FMIN)/(1 - (FMAX/FMIN)*exp(-x/omega));
end