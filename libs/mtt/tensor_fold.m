%% [A] = tensor_fold(Au,s)
%
function [A] = tensor_fold(Au,s)
  [n1] = s(1);
  [n2] = s(2);
  [n3] = s(3);
  A = zeros(n1,n2,n3);
  idx = 1;
  for i = 1:n3
    A(:,:,i) = Au(idx:idx+n1-1,:);
    idx = idx + n1;
  end
end
