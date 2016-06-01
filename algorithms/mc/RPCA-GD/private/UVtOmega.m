function y = UVtOmega(U,V,I,J,col); 
% The values of U*V' in Omega, as a vector
% Omega is given by I,J, col

y = zeros(length(I), 1);
for k = 1:length(col)-1
    j = J(col(k)+1);
    Xj = U * V(j,:)';
    idx = [col(k)+1:col(k+1)];
    y(idx) = Xj(I(idx));
end

% for j = 1:length(col)-1
%     Xj = U * V(j,:)';
%     idx = [col(j)+1:col(j+1)];
%     y(idx) = Xj(I(idx));
% end


