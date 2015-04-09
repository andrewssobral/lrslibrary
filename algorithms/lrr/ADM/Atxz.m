function Atz = Atxz(z);

global A Xg eta M;%A is the skinny SVD of Z_k, Xg is a copy of X, and M=X-E_{k+1}+Y/mu_k.

%compute Z'*z
temp1 = (A.U)'*z;
temp1 = (A.s).*temp1;
temp1 = (A.V)*temp1;

%compute X*z
Xz = Xg*z;

%compute X'*X*z = X'*(X*z)=X'*Xz
temp2 = Xg'*Xz;

%compute Z'*X'*X*z = Z'*(X'*X*z) = Z'*temp2
temp3 = (A.U)'*temp2;
temp3 = (A.s).*temp3;
temp3 = (A.V)*temp3;

%compute (X'*(X - X*Z - E - Y/mu))'*z = (M' - Z'*X')*(X*z)
% = M'*Xz - Z'*(X'*Xz) = M'*Xz - temp3
temp4 = M'*Xz -  temp3;

Atz = temp1 + temp4/eta;

