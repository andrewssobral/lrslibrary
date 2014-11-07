%% [A2_] = tensor_slices_lateral(A)
%
function [A2_] = tensor_slices_lateral(A)
  [n1, n2, n3] = size(A);
  for i = 1:n2
    a = A(:,i,:);
    A2_{i} = reshape(a,[n1 n3]);
  end
end
