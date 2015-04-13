function V = tens2vec(T,mode_row)
%TENS2VEC Vectorize a tensor.
%   V = tens2vec(T,mode_row) vectorizes a tensor T into a column vector V.
%   The rows of V are obtained by looping over the indices of T in the
%   order mode_row. E.g., if A and B are two matrices and T = cat(3,A,B),
%   then tens2vec(T,1:3) is the vector [A(:);B(:)].
%
%   V = tens2vec(T) vectorizes a tensor T, where mode_row is chosen as the
%   sequence 1:ndims(T).
%
%   See also vec2tens, tens2mat.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

if nargin < 2
    V = T(:);
else
    V = tens2mat(T,mode_row);
end
