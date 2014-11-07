function [lambda,x,flag,its,x0,trace,tracex] = sshopm(A,varargin)
%SSHOPM Shifted power method for finding a real eigenpair of a real tensor.
%
%   [LAMBDA,X]=SSHOPM(A) finds an eigenvalue (LAMBDA) and eigenvector (X)
%   for the real tensor A such that Ax^{m-1} = lambda*x.
%
%   [LAMBDA,X]=SSHOPM(A,parameter,value,...) can specify additional
%   parameters as follows: 
% 
%     'Shift'    : Shift in the eigenvalue calculation (Default: 0)
%     'MaxIts'   : Maximum power method iterations (Default: 1000)
%     'Start'    : Initial guess (Default: normal random vector)
%     'Tol'      : Tolerance on norm of change in |lambda| (Default: 1e-16)
%     'Concave'  : Treat the problem as concave rather than convex.
%                  (Default: true for negative shift; false otherwise.)
%     'Display'  : Display every n iterations (Default: -1 for no display)
%
%   [LAMBDA,X,FLAG]=SSHOPM(...) also returns a flag indicating convergence.
%
%      FLAG = 0  => Succesfully terminated 
%      FLAG = -1 => Norm(X) = 0
%      FLAG = -2 => Maximum iterations exceeded
%
%   [LAMBDA,X,FLAG,IT]=SSHOPM(...) also returns the number of iterations.
%
%   [LAMBDA,X,FLAG,IT,X0]=SSHOPM(...) also returns the intial guess.
%
%   [LAMBDA,X,FLAG,IT,X0,TRACE]=SSHOPM(...) also returns a trace of the
%   lambda values at each iteration.
%
%   REFERENCE: T. G. Kolda and J. R. Mayo, Shifted Power Method for
%   Computing Tensor Eigenpairs, SIAM Journal on Matrix Analysis and
%   Applications 32(4):1095-1124, October 2011 (doi:10.1137/100801482)    
%
%   See also SSHOPMC, TENSOR, SYMMETRIZE, ISSYMMETRIC.
%
%MATLAB Tensor Toolbox.
%Copyright 2012, Sandia Corporation.

% This is the MATLAB Tensor Toolbox by T. Kolda, B. Bader, and others.
% http://www.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2012) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in the file LICENSE.txt


%% Error checking on A
P = ndims(A);
N = size(A,1);

if ~issymmetric(A)
    error('Tensor must be symmetric.')
end

%% Check inputs
p = inputParser;
p.addParamValue('Shift', 0);
p.addParamValue('MaxIts', 1000, @(x) x > 0);
p.addParamValue('Start', [], @(x) isequal(size(x),[N 1]));
p.addParamValue('Tol', 1.0e-16);
p.addParamValue('Display', -1, @isscalar);
p.addParamValue('Concave', -1);
p.parse(varargin{:});

%% Copy inputs
maxits = p.Results.MaxIts;
x0 = p.Results.Start;
shift = p.Results.Shift;
tol = p.Results.Tol;
display = p.Results.Display;
concave = p.Results.Concave;

%% Check starting vector
if isempty(x0)
    x0 = 2*rand(N,1)-1;
end

if norm(x0) < eps
    error('Zero starting vector');
end

%% Check concavity
if concave == -1
    concave = (shift < 0);
end        

%% Execute power method
if (display >= 0)
    fprintf('TENSOR SHIFTED POWER METHOD: ');
    fprintf('Shift = %g, ', shift);
    fprintf('Concave = %d, ', concave);
    fprintf('\n');
    fprintf('----  --------- ----- --------\n');
    fprintf('Iter  Lambda    Diff  |newx-x|\n');
    fprintf('----  --------- ----- --------\n');
end

flag = -2;
x = x0 / norm(x0);
lambda = x'*ttsv(A,x,-1); 

trace = zeros(maxits,1);
trace(1) = lambda;
tracex = zeros(maxits,length(x));
tracex(1,:) = x0;

for its = 1:maxits
    
    newx = ttsv(A,x,-1) + shift * x;
    
    if (concave)
        newx = -newx;
    end
    
    nx = norm(newx);
    if nx < eps, 
        flag = -1; 
        break; 
    end
    newx = newx / nx;    
    
    newlambda = newx'* ttsv(A,newx,-1);    

    if norm(abs(newlambda-lambda)) < tol
        flag = 0;
    end
       
    if (display > 0) && ((flag == 0) || (mod(its,display) == 0))
        fprintf('%4d  ', its);
        % Lambda
        fprintf('%9.6f ', newlambda);
        d = newlambda-lambda;
        if (d ~= 0)
            if (d < 0), c = '-'; else c = '+'; end
            fprintf('%ce%+03d ', c, round(log10(abs(d))));
        else
            fprintf('      ');
        end
        % Change in X
        fprintf('%8.6e ', norm(newx-x));                
        % Line end
        fprintf('\n');
    end
    
    x = newx;
    lambda = newlambda;
    trace(its+1) = lambda;
    tracex(its+1,:) = x;
    
    if flag == 0
        break
    end
end

%% Check results
if (display >=0)
    switch(flag)
        case 0
            fprintf('Successful Convergence');
        case -1 
            fprintf('Converged to Zero Vector');
        case -2
            fprintf('Exceeded Maximum Iterations');
        otherwise
            fprintf('Unrecognized Exit Flag');
    end
    fprintf('\n');
end


%% ----------------------------------------------------
