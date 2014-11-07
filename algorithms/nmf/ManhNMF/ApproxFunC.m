function y = ApproxFunC(E,D,lambda)

% Closed-form Approximated Smoothing Function
% Written by Naiyang Guan (ny.guan@gmail.com)

X = abs(E)./D;
Y = X-lambda/2;
index = (X <= lambda);
Y(index) = X(index).^2/(lambda*2);
y = sum(sum(Y.*D));

return;