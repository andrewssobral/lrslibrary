%% [A1_] = tensor_slices_frontal(A)
%
function [A1_] = tensor_slices_frontal(A)
  [n1, n2, n3] = size(A);
  for i = 1:n3
    A1_{i} = A(:,:,i);
  end
end
