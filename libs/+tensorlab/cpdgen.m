function T = cpdgen(U)
%CPDGEN Generate full tensor given a polyadic decomposition.
%   T = cpdgen(U) computes the tensor T as the sum of R rank-one tensors
%   defined by the columns of the factor matrices U{n}.
%
%   See also btdgen, lmlragen, cpdres.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

T = reshape(U{1}*kr(U(end:-1:2)).',cellfun('size',U(:).',1));
