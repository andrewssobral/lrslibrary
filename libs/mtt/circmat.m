%% [M] = circmat(A)
%
function [Auc] = circmat(Au,s)
  [n1] = s(1);
  [n2] = s(2);
  [n3] = s(3);
  Auc = zeros(n1*n3,n2*n3);
  idx = 1;
  shi = n1;
  for k = 1:n3
    if k == 1
      Auc(:,idx:idx+n2-1) = Au(:,:);
    else
      Auc(:,idx:idx+n2-1) = circshift(Au,shi);
      shi = shi + n1;
    end
    idx = idx + n2;
  end
end
