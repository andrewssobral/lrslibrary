function out = drawFromIG(theta, chi)
    [m, n] = size(theta);
    chisq1 = randn(m, n).^2;
    out = theta + 0.5*theta./chi .* (theta.*chisq1 - sqrt(4*theta.*chi.*chisq1 + theta.^2.*chisq1.^2) );
    l = (rand(m, n) >= theta./(theta+out));
    out(l) = theta(l).^2 ./ out(l);
end