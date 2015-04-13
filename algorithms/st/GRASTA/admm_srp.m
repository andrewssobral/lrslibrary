function [  s, w, y ] = admm_srp( U, v, OPTS )
% Solve sparse residual pursuit (SRP) via ADMM {S.Boyd 2011}
%
% Disclaimer: This matlab function "admm_srp" is derived
% from the ADMM sovler for Least Absolute Deviations (LAD). We modified the
% original LAD function to adapt to our GRASTA framework. Interested readers
% can find other great matlab ADMM functions from S. Boyd's ADMM webpage.
%
% [ s, w, y, history ] = admm_srp( U, v, OPTS )
% 
% Solves the following problem via ADMM:
% 
%   minimize     ||s||_1
%   subject to   Uw + s- v = 0
%
% The solution is returned in the vectors s, w, y.
%
% y is the dual vector which will be used in GRASTA. 
%
% OPTS: structure for tuning the behavior of SRP
%   RHO   : the augmented Lagrangian parameter, default 1.8
%   MAX_ITER: max iteration for SRP problem
%   
%   We also provide MEX version for this function which improves the
%   performance at least 3 times.
%   
%
% Date: Aug. 10, 2012
%

%% Global constants and defaults

if isfield(OPTS,'RHO'),
    rho = OPTS.RHO;
else
    rho = 1.8;
end

if isfield(OPTS,'TOL'),
    TOL = OPTS.TOL;
else
    TOL = 1e-7;
end

if isfield(OPTS,'MAX_ITER'),
    MAX_ITER = OPTS.MAX_ITER;
else
    MAX_ITER = 50;
end

%% Data preprocessing

[m, n] = size(U);


w = zeros(n,1);
s = zeros(m,1);
y = zeros(m,1);
mu = 1.25/norm(v);

% precompute static variables for a-update (projection on to Ua=v-s)
P = (U'*U) \ (U');


%% ADMM solver

converge = false;
iter     = 0;

while ~converge && iter < MAX_ITER,
    iter = iter + 1;
    % w update
    w = P * (v -s - y/mu);
    
    % s update 
    Uw = U*w;
    s = shrinkage( v-Uw - y/mu, 1/mu); 

    % y update
    h = Uw + s - v;
    y = y + mu * h;

    mu = rho * mu; 
        
    % diagnostics, reporting, termination checks
    
    if (norm(h) < TOL),
        converge = true;    
    end
   
end

end


function y = shrinkage(a, kappa)
    y = max(0, a-kappa) - max(0, -a-kappa);
end