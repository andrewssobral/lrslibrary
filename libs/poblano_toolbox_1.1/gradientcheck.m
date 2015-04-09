function out = gradientcheck(FUN,x0,varargin)
%GRADIENTCHECK   Finite difference verification of analytic gradients.
%
%   OUT = GRADIENTCHECK(FUN,X0) computes the finite difference gradient
%   approximations to the gradient of FUN at X0. FUN is a handle for a
%   function that takes a single vector input and returns two arguments ---
%   the scalar function value and the vector-valued gradient. The output
%   contains the following informationfields:
%
%     OUT.G                 : analytic gradient (G)
%     OUT.GFD               : FD approximation of gradient (GFD)
%     OUT.MaxDiff           : maximum difference between G and GFD
%     OUT.MaxDiffInd        : index of maximum difference between G and GFD
%     OUT.NormGradientDiffs : 2-norm of G - GFD
%     OUT.GradientDiffs     : G - GFD
%     OUT.Params            : paramters used to compute FD approximations
%
%   OUT = GRADIENTCHECK(FUN,X0,'param,value) specifies a
%   parameters and its value. Available parameters are as follows:
%
%     DifferenceType - difference scheme to use {'forward'}
%       'forward'  : g_i = (f(x+h*e_i) - f(x))/h
%       'backward' : g_i = (f(x) - f(x-h*e_i))/h
%       'centered' : g_i = (f(x+h*e_i) - f(x-h*e_i))/(2h)
%
%     DifferenceStep - value of h in difference formulae {1e-8}
%
%
%   EXAMPLE
%
%   Suppose the function and gradient of the objective function are
%   specified in an mfile named mysin.m:
%
%     function [f,g]=example1(x,a)
%     if nargin < 2, a = 1; end
%     f = sum(sin(a*x));
%     g = a*cos(a*x);
%
%   We can call the gradient checker (using its default
%   parameters) using the command:
%
%     out = gradientcheck(@(x) example1(x,3), pi/4)
%
%   To change a parameter, we can specify a param/value input pair
%   as follows:
%
%     out = gradientcheck(@(x) example1(x,3), pi/4, 'DifferenceType', 'centered')
%
%   Alternatively, we can use a structure to define the parameters:
%  
%     params.DifferenceStep = 1e-6;
%     out = gradientcheck(@(x) example1(x,3), pi/4, params)
%
%MATLAB Poblano Toolbox.
%Copyright 2009-2012, Sandia Corporation.

%% Parse parameters

% Create parser
params = inputParser;

% Set parameters for this method
params.addParamValue('DifferenceType','forward', @(x) ismember(x,{'forward','backward','centered'}));
params.addParamValue('DifferenceStep',1e-8, @(x) x >= 1e-16);

% Parse input
params.parse(varargin{:});

%% Compute finite difference approximation of gradient
fdtype = params.Results.DifferenceType;
h = params.Results.DifferenceStep;

% setup variables
N = length(x0);
ei = zeros(size(x0)); 
g = zeros(N,1);
if ~(strcmp(fdtype,'backward')), fxp = g; end
if ~(strcmp(fdtype,'forward')), fxm = g; end

% function value at x0
[f,g] = feval(FUN,x0);

for i = 1:N
    ei(i) = h;
    if i>1
        ei(i-1)=0; 
    end
    if ~(strcmp(fdtype,'backward')), fxp(i) = feval(FUN,x0+ei); end
    if ~(strcmp(fdtype,'forward')), fxm(i) = feval(FUN,x0-ei); end
end

% Compute the FD approximations
switch fdtype
    case 'forward'
        gfd = (fxp-f)/h;      % foward difference
    case 'backward'
        gfd = (f-fxm)/h;      % backward difference
    case 'centered'
        gfd = (fxp-fxm)/2/h;  % central difference
end
g_gfd = g-gfd;

%% Output

% Analytic gradient and FD approximation
out.G = g;
out.GFD = gfd;

% Maximum difference (in magnitude)
[m,ind] = max(abs(g_gfd));
out.MaxDiff = g_gfd(ind);
out.MaxDiffInd = ind;

% Differences between gradient and finite difference gradient
out.NormGradientDiffs = norm(g_gfd);
out.GradientDiffs = g_gfd;

% Parameters used to compute the differences
out.Params = params.Results;
