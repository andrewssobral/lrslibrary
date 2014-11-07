function [res,K] = dotk(x,y,n)
%DOTK Dot product in K-fold precision.
%   dotk(x,y) computes the scalar product x'*y of the vectors x and y as in
%   K-fold precision, where K is chosen adaptively so that the result is
%   accurate up to machine precision. If x and y are N-D arrays, dotk(x,y)
%   operates along their first non-singleton dimension.
%
%   dotk(x,y,n) returns the dot product of x and y in the dimension n.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] T. Ogita, S. M. Rump, S. Oishi, "Accurate sum and dot product",
%       SIAM J. Sci. Comp., Vol. 26, No. 6, 2005, pp. 1955-1988.

% Permute x and y so that the n-th mode is the first one.
sx = size(x);
sy = size(y);
if any(sx ~= sy), error('dotk:xy','x and y should have equal size.'); end
if nargin < 3, n = find(sx > 1,1); end
if isempty(n) || size(x,n) == 1, res = conj(x).*y; K = 0; return; end

% Call recursively if either x or y is complex.
if ~isreal(x) && isreal(y)
    [resr,Kr] = dotk(real(x),y,n);
    [resi,Ki] = dotk(-imag(x),y,n);
    res = resr+resi*1i;
    K = Kr+Ki*1i;
    return;
elseif isreal(x) && ~isreal(y)
    [resr,Kr] = dotk(x,real(y),n);
    [resi,Ki] = dotk(x,imag(y),n);
    res = resr+resi*1i;
    K = Kr+Ki*1i;
    return;
elseif ~isreal(x) && ~isreal(y)
    x = cat(n,real(x),imag(x));
    [resr,Kr] = dotk(x,cat(n,real(y),imag(y)),n);
    [resi,Ki] = dotk(x,cat(n,imag(y),-real(y)),n);
    res = resr+resi*1i;
    K = Kr+Ki*1i;
    return;
end
if n > 1
    x = permute(x,[n 1:n-1 n+1:ndims(x)]);
    y = permute(y,[n 1:n-1 n+1:ndims(y)]);
end

% Initialize the algorithm.
res = zeros([2*sx(n) sx(1:n-1) sx(n+1:length(sx))]);
factor = 2^ceil(log2(2/eps(class(x)))/2)+1;

% Adaptive K-fold precision dot product.
[res(end,:),res(1,:)] = twoproduct(x(1,:),y(1,:),factor);
for i = 2:size(x,1)
    [tmp,res(i,:)] = twoproduct(x(i,:),y(i,:),factor);
    [res(end,:),res(i+size(x,1)-1,:)] = twosum(res(end,:),tmp);
end
[res,K] = sumk(res);
if n > 1, res = permute(res,[2:n 1 n+1:length(sx)]); end
if isreal(K), K = K+1; else K = K+1+1*1i; end

end

function [x,y] = twoproduct(a,b,factor)

tmp = factor*a;
a1 = tmp-(tmp-a);
a2 = a-a1;
tmp = factor*b;
b1 = tmp-(tmp-b);
b2 = b-b1;

x = a.*b;
y = a2.*b2-(((x-a1.*b1)-a2.*b1)-a1.*b2);

end

function [x,y] = twosum(a,b)

x = a+b;
tmp = x-a;
y = (a-(x-tmp))+(b-tmp);

end
