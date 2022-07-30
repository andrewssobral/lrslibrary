function val = tensor_1norm( X )

X = double(X);
val = sum( abs(X(:)) );

end