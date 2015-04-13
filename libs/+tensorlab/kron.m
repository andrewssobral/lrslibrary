function X = kron(A,B)
%KRON Kronecker product.
%   kron(A,B) returns the Kronecker product of two matrices A and B, of 
%   dimensions I-by-J and K-by-L respectively. The result is an I*K-by-J*L
%   block matrix in which the (i,j)-th block is defined as A(i,j)*B.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

[I,J] = size(A);
[K,L] = size(B);

if ~issparse(A) && ~issparse(B)
    
    % Both matrices are dense.
    A = reshape(A,[1 I 1 J]);
    B = reshape(B,[K 1 L 1]);
    X = reshape(bsxfun(@times,A,B),[I*K J*L]);
    
else
    
    % One of the matrices is sparse.
    [ia,ja,sa] = find(A);
    [ib,jb,sb] = find(B);
    ix = bsxfun(@plus,K*(ia(:)-1).',ib(:));
    jx = bsxfun(@plus,L*(ja(:)-1).',jb(:));
    if islogical(sa) && islogical(sb)
        X = sparse(ix,jx,bsxfun(@and,sb(:),sa(:).'),I*K,J*L);
    else
        X = sparse(ix,jx,double(sb(:))*double(sa(:).'),I*K,J*L);
    end
    
end
