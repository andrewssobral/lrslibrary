function dFdz = deriv(F,z0,Fz0,type)
%DERIV Approximate gradient and Jacobian.
%   For real x0, deriv(F,x0) computes the gradient dF(x0)/dx or Jacobian
%   dF(x0)/dx^T of a real-valued function F(x) in the real variables x,
%   assuming F is analytic in x when x is complex. If this assumption is
%   not satisfied (e.g., when f is function of x' instead of x.'), use the
%   finite difference methods below. The input variables x may be a scalar,
%   vector, matrix, tensor or even a (nested) cell array of tensors.
%
%   For complex z0, deriv(F,z0) and deriv(F,z0,Fz0) are equivalent to
%   deriv(F,z0,Fz0,'gradient') if F(z0) is a real scalar, and
%   deriv(F,z0,Fz0,'Jacobian') otherwise. The input variables z may be a
%   scalar, vector, matrix, tensor or even a cell array of tensors.
%
%   deriv(f,z0,fz0,'gradient') approximates the real gradient df(z0)/dx or
%   scaled conjugate cogradient 2*conj(df(z0)/dz) of a real-valued function
%   f in the real variables x or complex variables z, respectively. The
%   scaled conjugate cogradient is defined as twice the complex conjugate
%   of the partial derivative of f w.r.t. the complex variables z, while
%   treating conj(z) as constant. The real gradient is returned if z0 is
%   real, else the scaled conjugate cogradient is returned. Note that the
%   former is equal to the latter if both z0 and f(z) are real. The input
%   fz0 should supply f(z0). If fz0 is [], it is computed as f(z0).
%
%   deriv(F,z0,Fz0,'Jacobian') approximates the Jacobian dF(z0)/dz^T in the
%   real or complex variables z. If F(z) is nonanalytic in z and depends on
%   conj(z), use 'Jacobian-C' instead. The input Fz0 should supply F(z0).
%   If Fz0 is [], it is computed as F(z0).
%
%   deriv(F,z0,Fz0,'Jacobian-C') approximates the complex Jacobian
%   [dF(z0)/dz^T dF(z0)/dconj(z)^T] of a function F in the complex
%   variables z. The partial derivatives in the complex Jacobian are w.r.t.
%   the complex variables z and conj(z), while treating conj(z) and z as
%   constant, respectively. The input Fz0 should supply F(z0). If Fz0 is
%   [], it is computed as F(z0).
%
%   See also gradient, diff.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Unconstrained
%       optimization of real functions in complex variables", SIAM J. Opt.,
%       Vol. 22, No. 3, 2012, pp. 879-898.

% If possible, use the i-trick.
dim = structure(z0);
z0 = serialize(z0);
F = @(z)serialize(F(deserialize(z,dim)));
if nargin < 3 && all(isreal(z0))
    % Assume F is analytic and both F and z0 are real-valued; return dF/dx.
    e = 2^floor(log2(nthroot(realmin(class(z0)),3)));
    p = zeros(size(z0));
    p(1) = e;
    dFdz = imag(F(z0+p*1i))/e;
    dFdz = [dFdz zeros(numel(dFdz),length(z0)-1)];
    for n = 2:length(z0)
        p(n-1) = 0; p(n) = e;
        dFdz(:,n) = imag(F(z0+p*1i))/e;
    end
    if size(dFdz,1) == 1
        dFdz = deserialize(dFdz(:),dim);
    end
    return;
elseif nargin < 3 || isempty(Fz0)
    Fz0 = F(z0);
end
Fz0 = serialize(Fz0);

% Define type of derivative.
if nargin < 4
    if numel(Fz0) == 1 && isreal(Fz0), type = 'gradient';
    else type = 'Jacobian'; end
end
dFdzIsGrad = any(strfind(type,'gradient'));
dFdzIsRealGrad = dFdzIsGrad && all(isreal(z0));
dFdzIsAnaJac = ~dFdzIsGrad && ~any(strfind(type,'-C'));

% Else use finite (complex) derivatives.
e = sqrt(eps(class(z0)));
p = zeros(size(z0));
p(1) = e*max(1,abs(real(z0(1))));
dFdx = (F(z0+p)-Fz0)/p(1);
if ~dFdzIsRealGrad || ~dFdzIsAnaJac
    p(1) = e*max(1,abs(imag(z0(1))))*1i;
    dFdy = (F(z0+p)-Fz0)/imag(p(1));
end
switch strtok(type,'-')
    case 'gradient'
        if dFdzIsRealGrad, dFdz = dFdx;
        else dFdz = real(dFdx)+real(dFdy)*1i; end
        dFdz = [dFdz zeros(numel(dFdz),length(z0)-1)];
    case 'Jacobian'
        if dFdzIsAnaJac
            dFdz = dFdx;
            dFdz = [dFdz zeros(numel(dFdz),length(z0)-1)];
        else
            dFdz = 0.5*cat(3,dFdx-dFdy*1i,dFdx+dFdy*1i);
            dFdz = cat(2,dFdz,zeros(size(dFdz,1),length(z0)-1,2));
        end
end
for n = 2:length(z0)
    p(n-1) = 0; p(n) = e*max(1,abs(real(z0(n))));
    dFdx = (F(z0+p)-Fz0)/p(n);
    if ~dFdzIsRealGrad || ~dFdzIsAnaJac
        p(n) = e*max(1,abs(imag(z0(n))))*1i;
        dFdy = (F(z0+p)-Fz0)/imag(p(n));
    end
    switch strtok(type,'-')
        case 'gradient'
            if dFdzIsRealGrad, dFdz(:,n) = dFdx;
            else dFdz(:,n) = real(dFdx)+real(dFdy)*1i; end
        case 'Jacobian'
            if dFdzIsAnaJac, dFdz(:,n) = dFdx;
            else dFdz(:,n,:) = 0.5*cat(3,dFdx-dFdy*1i,dFdx+dFdy*1i); end
    end
end
dFdz = dFdz(:,:);
if size(dFdz,1) == 1 && dFdzIsGrad
    dFdz = deserialize(dFdz(:),dim);
end

end

function [z,offset] = deserialize(z,dim,offset)
    if iscell(dim)
        v = z;
        z = cell(size(dim));
        if nargin < 3, offset = 0; end
        for i = 1:numel(z)
            if iscell(dim{i})
                [z{i},offset] = deserialize(v,dim{i},offset);
            else
                n = prod(dim{i}(:));
                z{i} = reshape(v(offset+(1:n)),dim{i});
                offset = offset+n;
            end
        end
    elseif ~isempty(dim)
        z = reshape(z,dim);
    end
end

function z = serialize(z)
    if iscell(z)
        for i = find(cellfun(@iscell,z(:).'))
            z{i} = serialize(z{i});
        end
        s = cellfun(@numel,z(:)); o = [0; cumsum(s)];
        c = z; z = zeros(o(end),1);
        for i = 1:length(s), z(o(i)+(1:s(i))) = c{i}(:); end
    else
        z = z(:);
    end
end

function dim = structure(z)
    if iscell(z)
        dim = cellfun(@size,z,'UniformOutput',false);
        for i = find(cellfun(@iscell,z(:).'))
            dim{i} = structure(z{i});
        end
    else
        dim = size(z);
        if numel(z) == dim(1), dim = []; end
    end
end
