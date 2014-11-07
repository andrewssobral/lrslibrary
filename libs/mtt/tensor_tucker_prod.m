%% [B] = tensor_tucker_prod(A,M)
%
function [B] = tensor_tucker_prod(A,M)
  B = A;
  for k = 1:length(A.size)
    B = ttm(B, M{k}, k);
  end
end
