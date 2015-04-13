function X = kr(U,varargin)
%KR Khatri-Rao product.
%   kr(A,B) returns the Khatri-Rao product of two matrices A and B, of 
%   dimensions I-by-K and J-by-K respectively. The result is an I*J-by-K
%   matrix formed by the matching columnwise Kronecker products, i.e.,
%   the k-th column of the Khatri-Rao product is defined as
%   kron(A(:,k),B(:,k)).
%
%   kr(A,B,C,...) and kr({A B C ...}) compute a string of Khatri-Rao 
%   products A x B x C x ..., where x denotes the Khatri-Rao product.
%
%   See also kron.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

if ~iscell(U), U = [{U} varargin]; end
[J,K] = size(U{end});
if any(cellfun('size',U,2) ~= K)
    error('kr:U','Input matrices should have the same number of columns.');
end

X = reshape(U{end},[J 1 K]);
for n = length(U)-1:-1:1
    I = size(U{n},1);
    A = reshape(U{n},[1 I K]);
    X = reshape(bsxfun(@times,A,X),[I*J 1 K]);
    J = I*J;
end
X = reshape(X,[size(X,1) K]);
