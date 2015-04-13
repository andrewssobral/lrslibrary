function [res,K] = sumk(p,n)
%SUMK Summation in K-fold precision.
%   sumk(p) sums a vector p as in K-fold precision, where K is chosen
%   adaptively so that the result is accurate up to machine precision. If p
%   is an N-D array, sumk(p) operates along the first non-singleton
%   dimension.
%
%   sumk(p,n) sums along the dimension n.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] T. Ogita, S. M. Rump, S. Oishi, "Accurate sum and dot product",
%       SIAM J. Sci. Comp., Vol. 26, No. 6, 2005, pp. 1955-1988.

% Permute p so that the n-th mode is the first one.
sp = size(p);
if nargin < 2, n = find(sp > 1,1); end
if isempty(n) || size(p,n) == 1, res = p; K = 0; return; end

% Call recursively if p is complex.
if ~isreal(p)
    [resr,Kr] = sumk(real(p),n);
    [resi,Ki] = sumk(imag(p),n);
    res = resr+resi*1i;
    K = Kr+Ki*1i;
    return;
end
if n > 1, p = permute(p,[n 1:n-1 n+1:ndims(p)]); end

% Initialize the algorithm.
res = inf([sp(1:n-1) 1 sp(n+1:length(sp))]);
K = 0;
idx = 1:numel(res);
stop = false(1,numel(res));

% Adaptive K-fold precision sum.
while ~all(stop);
    K = K+1;
    for i = idx(~stop)
        [p(:,i),resi] = vecsum(p(:,i));
        stop(i) = resi >= res(i) || resi < eps(p(end,i));
        res(i) = resi;
    end
end
res = reshape(sum(p(1:end-1,:),1)+p(end,:),[1 sp([1:n-1 n+1:length(sp)])]);
if n > 1, res = permute(res,[2:n 1 n+1:length(sp)]); end

end

function [p,ce] = vecsum(p)

ce = 0;
for i = 2:length(p)
    a = p(i);
    b = p(i-1);
    p(i) = a+b;
    tmp = p(i)-a;
    p(i-1) = (a-(p(i)-tmp))+(b-tmp);
    ce = max(ce,abs(p(i-1)));
end
ce = (length(p)-1)*ce;

end
