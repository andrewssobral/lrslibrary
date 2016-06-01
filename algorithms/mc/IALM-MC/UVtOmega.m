function y = UVtOmega(U,V,I,J,col); 
%Computes the restriction of UV' onto Omega, which is defined by the list (I,J)

y = zeros(length(I), 1);
for k = 1:length(col)-1
    j = J(col(k)+1);
    Xj = U * V(j,:)';
    idx = [col(k)+1:col(k+1)];
    y(idx) = Xj(I(idx));
end

