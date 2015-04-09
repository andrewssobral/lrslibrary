function [f,g]=example1(x,a)
%EXAMPLE1   Example function for Poblano Toolbox routines.
%
%  [F,G] = EXAMPLE1(X) computes F = sin(X(N)) + ... sin(X(N)) and G, the 
%  gradient of F at X, where N is the number of variables in X.
%
%  [F,G] = EXAMPLE1(X,A) computes F = sin(A*X(1)) + ... + sin(A*X(N)) and G, 
%  the gradient of F at X.
%
%MATLAB Poblano Toolbox.
%Copyright 2009-2012, Sandia Corporation.

if (nargin < 2), a = 1; end
f = sum(sin(a*x));
g = a*cos(a*x);

