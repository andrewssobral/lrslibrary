function tnorm = tensor_array_norm( X )
% X is a cell array of tensors

tnorm = 0;
for i = 1:length(X)
    tnorm = tnorm + norm(X{i})^2;
end
tnorm = sqrt(tnorm);

end