%%
% The lansvd function can be downloaded here
% http://soi.stanford.edu/~rmunk/PROPACK/
function [A,E,tol_iter] = alm(D)

[n,p] = size(D);
kappa = 1.1;%can be tuned
tau = 0.61;%can be tuned
lambda = tau*kappa;
eta = (1-tau)*kappa;
mu = 30/norm(sign(D)); 
rho = 1.1; 
alpha = 1;
beta = .2;

Y = zeros(size(D));
E = zeros(size(D));
A = D;

tol_inner1 = 1e-4;
tol_inner2 = 1e-6;
tol_out = 1e-7;
err_out = 1;
MaxIter_out = 500;
MaxIter_inner1 = 1;
MaxIter_inner2 = 20;

iter_out = 0;
sv = 10;
tol_iter = 0;

while iter_out < MaxIter_out && err_out > tol_out
    iter_out = iter_out + 1; disp(iter_out);

    Ak = A; Ek = E;
    iter_inner1 = 0;
    err_inner1 = 1;

    while iter_inner1 < MaxIter_inner1 && err_inner1 > tol_inner1
        iter_inner1 = iter_inner1 + 1;

        G = D - Ek + Y/mu;
        Akk = G;
        Ahk = zeros(size(Akk));

        err_inner2 = 1;
        iter_inner2 = 0;

        while iter_inner2 < MaxIter_inner2 && err_inner2 > tol_inner2
            iter_inner2 = iter_inner2 + 1;
            
            %if choosvd(n, sv) == 1
            %    [U,Si,V] = lansvd(Akk, sv, 'L');
            %else
                %[U,Si,V] = svd(Akk, 'econ');
                [U,Si,V] = svdecon(Akk); % fastest
            %end
            
            diagS = diag(Si);
            diagS = diagS(1:sv);
            svn = length(find(diagS > beta));
            svp = svn;
            ratio = diagS(1:end-1)./diagS(2:end);
            [max_ratio, max_idx] = max(ratio);
            if max_ratio > 2
                svp = min(svn, max_idx);
            end
            if svp < sv
                sv = min(svp + 1, n);
            else
                sv = min(svp + 10, n);
            end

            Ahk = U(:, 1:svp) * diag(diagS(1:svp) - beta) * V(:, 1:svp)';

            B = 2*Ahk - Akk + mu*beta*G;
            %ns = norms(B);
            ns = norm(B);
            B = bsxfun(@times,B/(1 + mu*beta),max(0,1 - beta*eta./ns));
            Akk = Akk + alpha*(B - Ahk);
            err_inner2 = alpha*norm(B-Ahk,'fro');

            tol_iter = tol_iter + 1;
        end

        G = D - Ahk + Y/mu;
        %ns = norms(G);
        ns = norm(G);
        Ep = bsxfun(@times, G, max(0,1-lambda/mu./ns));

        err_inner1 = max(norm(Ek-Ep,'fro'),norm(Ak-Ahk,'fro'));
        Ek = Ep;
        Ak = Ahk;
    end

    A = Ak; E = Ek;
    err_out = norm(D-A-E,'fro')/norm(D,'fro');

    Y = Y + mu*(D - A - E);
    mu = rho*mu;
end

