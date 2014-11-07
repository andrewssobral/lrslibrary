%% [At] = tensor_transpose(A)
%  p-transpose, p<[2 1 3]>
function [At] = tensor_transpose(A)
  At = permute(A, [2 1 3]);
end
