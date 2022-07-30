function D = tensor_array_diff( X1, X2 )

D = cell( 1, length(X1) );
for i = 1:length(X1)
    D{i} = X1{i} - X2{i};
end
end