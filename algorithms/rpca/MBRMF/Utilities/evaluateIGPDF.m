function res = evaluateIGPDF(mu, lambda, x)
    res = 0.5 * log(lambda ./ (2 * pi * x .^ 3)) + (-lambda .* (x - mu) .^ 2 ./ (2 * (mu .^ 2) .* x));
end