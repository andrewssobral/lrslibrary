%*******************************************************************
% Euclidean distance between two matrixs' column vector
%*******************************************************************
function [D] = EucliDist2(A, B)

if size(A, 1) == size(B, 1)
    D = A'.^2*ones(size(B))+ones(size(A'))*(B).^2-2*A'*B;
else
    disp('incorrect input matrix, the first dimension not matched');
    D = 0;
end