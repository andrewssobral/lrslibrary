function X = normalize(X)
    [sizeD1, sizeD2] = size(X);
    X = X - ones(sizeD1, 1) * mean(X);
    DTD = X' * X;
    invTrX = ones(sizeD2, 1) ./ sqrt(diag(DTD));
    mul = ones(sizeD1, 1) * invTrX';
    X =X .* mul;
end