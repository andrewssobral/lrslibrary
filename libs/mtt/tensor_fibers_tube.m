%% [B3_] = tensor_fibers_tube(A)
%
function [B3_] = tensor_fibers_tube(A)
  [n1, n2, n3] = size(A);
  for i = 1:n1
    for j = 1:n2
      B3_{i,j} = A(i,j,:);
    end
  end
end
