%% [B1_] = tensor_fibers_column(A)
%
function [B1_] = tensor_fibers_column(A)
  [n1, n2, n3] = size(A);
  for j = 1:n2
    for k = 1:n3
      B1_{j,k} = A(:,j,k);
    end
  end
end
