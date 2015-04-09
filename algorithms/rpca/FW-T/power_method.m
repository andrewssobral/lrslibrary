function [ u s v ] = power_method( A, k )

% by power method
[m n] = size(A);

if k<=10
    if m<=n
        q = randn(m,1);
        for i = 1:k
            q = A*(A'*q);
            q = q/norm(q,2);
        end
        u = q; v = A'*q;
        s = norm(v,2);
        v = v/s;
    else
        q = randn(n,1);
        for i = 1:k
            q = A'*(A*q);
            q = q/norm(q,2);
        end
        v = q; u = A*v;
        s = norm(u,2);
        u = u/s;
        
    end
end

if k>10
    if m<=n
        H = A*A'; q = randn(m,1);
        for i = 1:k
            q = H*q;
            q = q/norm(q,2);
        end
        s = sqrt(q'*H*q);
        u = q;
        v = A'*u/s;
    else
        H = A'*A; q = randn(n,1);
        for i = 1:k
            q = H*q;
            q = q/norm(q,2);
        end
        s = sqrt(q'*H*q);
        v = q;
        u = A*v/s;
    end
end
end

