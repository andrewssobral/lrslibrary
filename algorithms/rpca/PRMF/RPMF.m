function [U V] = RPMF(Y, rk, lambdaU, lambdaV, tol)
    %% Initialize essencial variables
    maxIter = 40;
    [m n] = size(Y);
    U = randn(m, rk);
    V = randn(rk, n);
    lambda = 1;
    eps = 1e-3;
    r = abs(Y - U * V);
    r = (r < eps) * eps + (r > eps) .* r;
    r = (lambda) ./ r;
    c = 0;
    IS = sparse(eye(rk));
    %% Main Process
    %tic;
    while true
        c = c + 1;
        disp(c);
        oldR = r;
       %% Update V
            parfor i = 1 : n
                T = (U' * diag(sparse(r(:, i))));
                    V(:, i) = (T * U + (lambdaV)* IS) \ (T * Y(:, i));
                r(:, i) = abs(Y(:, i) - U * V(:, i));
                r(:, i) = (r(:, i) < eps) * eps + (r(:, i) > eps) .* r(:, i);
                r(:, i) = (lambda) ./ r(:, i);
            end
       %% Update U
        parfor i = 1 : m
            T = (V * diag(sparse(r(i, :))));
            U(i, :) = (T * V' + lambdaU * IS) \ (T * Y(i, :)');
            r(i, :) = abs(Y(i, :)' - V' * U(i, :)');
            r(i, :) = (r(i, :) < eps) * eps + (r(i, :) > eps) .* r(i, :);
            r(i, :) = (lambda) ./ r(i, :);
        end
        if sum(abs(r(:) - oldR(:))) / sum(oldR(:)) < tol && c ~= 1 || c > maxIter
            break;
        end
    end
    %toc
end