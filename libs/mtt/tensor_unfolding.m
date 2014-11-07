%% [An] = tensor_unfolding(A)
% matricization
function [An] = tensor_unfolding(A,mode)
  An = double(tenmat(tensor(A),mode));
end
