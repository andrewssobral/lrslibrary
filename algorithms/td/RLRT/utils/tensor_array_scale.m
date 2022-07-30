function X = tensor_array_scale( X, a )

for i = 1:length(X)
    X{i} = X{i}*a;
end

end