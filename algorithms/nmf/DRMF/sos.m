function [result]=sos(A, dim)
%[result] = sos(A, dim)
% return the sum of squares
% if dim = 0, then return the total sos
% otherwise, as in sum or mean
% author: Liang Xiong (lxiong@cs.cmu.edu)

if nargin < 2
    dim=0;
end

if dim == 0
    result = sum(A(:).^2);
else
    result = sum(A.^2, dim);
end
