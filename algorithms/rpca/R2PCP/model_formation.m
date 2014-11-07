vm.f=@(A,B).5*norm(A(:)+B(:)-Z(:))^2;
vm.g=@(A,B)A+B-Z;
