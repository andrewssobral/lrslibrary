function [W,H] = nndSVD(A, k)
nSample = size(A, 1);
nFeature = size(A, 2);
k = k + 1;
% --------------------------------------------
[U, S, V] = svds(A, k);

% --------------------------------------------
W = zeros(nSample, k);
H = zeros(k, nFeature);
W(:,1) = sqrt(S(1,1))*U(:,1);
H(1,:) = sqrt(S(1,1))*V(:,1)';
for i = 2:k
    x = U(:,i);
    y = V(:,i);
    xp = pos(x); xn = neg(x);
    yp = pos(y); yn = neg(y);
    xpnorm = norm(xp); ypnorm = norm(yp); mp = xpnorm*ypnorm;
    xnnorm = norm(xn); ynnorm = norm(yn); mn = xnnorm*ynnorm;
    if mp > mn
        u = xp/xpnorm;
        v = yp/ypnorm;
        sigma = mp;
    else
        u = xn/xnnorm;
        v = yn/ynnorm;
        sigma = mn;
    end
    W(:,i) = sqrt(S(i,i)*sigma)*u;
    H(i,:) = sqrt(S(i,i)*sigma)*v';
end
W = W(:, 2:k);
H = H(2:k, :);
end

function [xp] = pos(x)
xp = (x >= 0).*x;
end

function [xn] = neg(x)
xn = (x < 0).*(-x);
end
