function f = frob(T,squared)
%FROB Frobenius norm.
%   frob(T) returns the Frobenius norm of the tensor T. If T is an
%   incomplete tensor, the Frobenius norm of the known elements is
%   returned.
%
%   frob(T,'squared') returns the squared Frobenius norm of the tensor T.
%
%   See also norm.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

if nargin < 2, squared = false; end
if isstruct(T)
    if ischar(squared), f = T.val'*T.val;
    else f = norm(T.val); end
else
    if ischar(squared), f = T(:)'*T(:);
    else f = norm(T(:)); end
end
