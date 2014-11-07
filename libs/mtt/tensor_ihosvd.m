function [T_hat] = tensor_ihosvd(core, U)
  T_hat = core;
  
  for i = 1:ndims(core)
    %[m n] = size(U{i});
    %s(i) = m;
    %T_hat = tensor(folding(U{i}*unfolding(double(T_hat),i), i, s));
    T_hat = ttm(T_hat,U{i},i);
  end
end
