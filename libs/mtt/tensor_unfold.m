%% [Au] = tensor_unfold(A)
%
function [Au] = tensor_unfold(A)
  [n1, n2, n3] = size(A);
  Au = zeros(n1*n3,n2);
  idx = 1;
  for i = 1:n3
    A_{i} = A(:,:,i);
    Au(idx:idx+n1-1,:) = A_{i};
    idx = idx + n1;
  end
end
