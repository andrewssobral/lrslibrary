function X = tensor_array_add( X1, X2 )
N = length(X1);
X = cell(1,N);
for i = 1:N
    X{i} = X1{i} + X2{i};
end
end