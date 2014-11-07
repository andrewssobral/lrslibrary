function [i,j] = argmin2(x)
%ARGMIN2  Index of minimum element of matrix.
% [i,j] = ARGMIN2(x) returns indices (i,j) such that x(i,j) == min(x(:)).
%
% See also ARGMIN, ARGMAX2.

[colmin,i] = min(x);
[ignore,j] = min(colmin);
i = i(j);
