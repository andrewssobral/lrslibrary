%%**************************************************
% compute A*z for 
% A = A_U * A_V' + Sparse_Z;
% 
% Atz = matvec(z,A_U,A_V,Sparse_Z); 
%
%%**************************************************
%%
  function Atz = Atxz(z); 

  global A Sparse_Z
  
  Atz = A.V * (A.U' * z) + (z' * Sparse_Z)';