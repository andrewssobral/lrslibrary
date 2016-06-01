%%**************************************************
% compute A*z for 
% A = A_U * A_V' + Sparse_Z;
% 
% Az = matvec(z,A_U,A_V,Sparse_Z); 
%
%%**************************************************
%%
  function Az = Axz(z); 

  global A Sparse_Z
  
  Az = A.U * (A.V' * z) + Sparse_Z * z;