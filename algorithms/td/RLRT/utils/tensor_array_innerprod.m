function prod = tensor_array_innerprod( X1, X2 )

prod = 0;
N = length(X1);
for i = 1:N
    prod = prod + innerprod(X1{i},X2{i});
end

end