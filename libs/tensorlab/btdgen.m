function T = btdgen(U)
%BTDGEN Generate full tensor given a BTD.
%   T = btdgen(U) computes the tensor T as the sum of the R block terms
%   U{r} which are computed as tmprod(U{r}{N+1},U{r}(1:N),1:N), where N is
%   the order of the tensor T.
%
%   See also cpdgen, lmlragen, btdres.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

N = length(U{1})-1;
size_tens = cellfun('size',U{1}(1:N),1);
T = reshape(U{1}{1}*mtkronprod(U{1}{end},U{1}(1:N),1,'H'),size_tens);
for r = 2:length(U)
    T = T+reshape(U{r}{1}*mtkronprod(U{r}{end},U{r}(1:N),1,'H'),size_tens);
end
