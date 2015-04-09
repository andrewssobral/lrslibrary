function Az = Axz(z);

global A Xg eta M;%A is the skinny SVD of Z_k, Xg is a copy of X, and M=X-E_{k+1}+Y/mu_k.

%compute Z*z
temp1 = (A.V)'*z;
temp1 = (A.s).*temp1;
temp1 = A.U*temp1;

%compute (X - X*Z - E - Y/mu)*z = (M - X*Z)z = M*z - X*(Z*z)
temp2 = M*z - Xg*temp1;

%compute X'*(X - X*Z - E - Y/mu)*z
temp3 = Xg'*temp2;

Az = temp1 + temp3/eta;

