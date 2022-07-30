% switch between the lansvd and standard svd, fast implementation
function [U S V] = smartSVD(A, k)

[m n] = size(A);
dim = min(m,n);
       
if 10*m < n
    AAT = A*A';
    
    if k==dim
        [U, S, V] = svd(AAT); 
    else
        [U, S, V] = lansvd(AAT, k, 'L');
    end
    S = diag(S) .^ 0.5;

    tol = max(size(A)) * eps(max(S));
    R = sum(S > tol);
    U = U(:,1:R);
    S = S(1:R);

    V = A'*U*diag(1./S);
    S = diag(S);
    return;
end
    
if m > 10*n
    ATA = A'*A;
    
    if k==dim
        [U, S, V] = svd(ATA); 
    else
        [U, S, V] = lansvd(ATA, k, 'L');
    end
    S = diag(S) .^ 0.5;

    tol = max(size(A)) * eps(max(S));
    R = sum(S > tol);
    V = V(:,1:R);
    S = S(1:R);

    U = A*V*diag(1./S);
    S = diag(S);
    return;
end

if k==dim
    [U, S, V] = svd(A);
else
    [U, S, V] = lansvd(A, k, 'L');
end
