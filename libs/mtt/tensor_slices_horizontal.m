%% [A3_] = tensor_slices_horizontal(A)
%
function [A3_] = tensor_slices_horizontal(A)
  [n1, n2, n3] = size(A);
  for i = 1:n1
    a = A(i,:,:);
    A3_{i} = reshape(a,[n2 n3]);
  end
end
