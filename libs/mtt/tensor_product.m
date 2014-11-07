function [C] = tensor_product(A,B)
  Au = tensor_unfold(A);
  Bu = tensor_unfold(B);
  
  Auc = circmat(Au,size(A));
  
  AucBu = Auc*Bu;
  
  C = tensor_fold(AucBu,size(A));
end

