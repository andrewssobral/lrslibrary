function T = lmlragen(U,S)
%LMLRAGEN Generate full tensor given a core tensor and factor matrices.
%   T = lmlragen(U,S) computes the tensor T as the mode-n tensor-matrix
%   product of the factor matrices U{n} with the core tensor S.
%
%   See also btdgen, cpdgen, lmlrares.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

T = btdgen({[U(:).',S]});
