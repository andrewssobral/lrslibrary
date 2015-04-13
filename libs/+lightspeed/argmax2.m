function [i,j] = argmax2(x)
%ARGMAX2  Index of maximum element of matrix.
% [i,j] = ARGMAX2(x) returns indices (i,j) such that x(i,j) == max(x(:)).
%
% See also ARGMAX.

[colmax,i] = max(x);
[ignore,j] = max(colmax);
i = i(j);
