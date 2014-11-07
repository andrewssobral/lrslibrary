%% [B2_] = tensor_fibers_row(A)
%
function [B2_] = tensor_fibers_row(A)
  [n1, n2, n3] = size(A);
  for i = 1:n1
    for k = 1:n3
      B2_{i,k} = A(i,:,k);
    end
  end
end
