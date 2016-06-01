% function out = getoptS_mex(M,X,Y,A,I,J)
% compute min_S |P_M(XSY')-A|_F
% M: sparse mask matrix (m x n) (values will be changed due to the mex trick)
% A: non-zero entris (nnz x 1)
% I,J: row, col indices (nnz x 1)
% Adapted from Keshavan et al.
function out = getoptS_mex(M,X,Y,A,I,J,lam)

[m r] = size(X);
[n r] = size(Y);
l = length(I);

C = X'*multsparsefull(A,Y,M);
C = C(:);

for i = 1:r
  for j = 1:r
    ind = (j-1)*r + i;
    temp = X'*multsparsefull(maskmult(X(:,i)',Y(:,j)',I,J),Y,M);
    SS(:,ind) = temp(:) ;
    if lam>0.0
      SS(ind,ind) = SS(ind,ind) + lam;
    end
  end
end
S = SS\C ;
out = reshape(S,r,r);

