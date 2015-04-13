function [xy,res,dxy] = polysol2(p,q,options)
%POLYSOL2 Solve a system of two bivariate polynomials.
%   [xy,res,dxy] = polysol2(p,q) computes isolated solutions xy of the
%   polynomial system p = q = 0, where p and q respresent polynomials in
%   two real variables (x,y) or in one complex variable z (and conj(z)).
%   The matrices p and q represent one of the following polynomial types:
%
%      [1 y  y^2  ... y^dy] *p*[1;x;x^2;...;x^dx]     (analytic  bivariate)
%      [1 z' z'^2 ... z'^d1]*p*[1;z;z^2;...;z^d2] (polyanalytic univariate)
%
%   where p(i,j) is the coefficient of the term x^j*y^i and z^j*conj(z)^i,
%   respectively. Each row of xy(i,:) is a solution of p = q = 0, res(i,:)
%   is equal to the normalized residual
%
%      [p(xy(i,:))/|p|(|xy(i,:)|) q(xy(i,:))/|q|(|xy(i,:)|)],
%
%   and dxy(:,:,i) is the real or complex Jacobian [dp/dx dp/dy;
%   dq/dx dq/dy] or [dp/dxy dp/d(conj(z)); dq/dxy dq/d(conj(z))],
%   respectively, at xy(i,:). If p and q are bivariate polynomials in x and
%   y, only real solutions of p(x,y) = q(x,y) = 0 are computed.
%
%   polysol2(p,options) may be used to set the following options:
%
%      options.MaxIter = 8      - The maximum number of iterations in the
%                                 Newton-Raphson refinement.
%      options.Plot = false     - If true, a plot is made of the contour
%                                 lines of level 0 for p and q.
%      options.PlotPoints =     - The number of grid points per axis to use
%      1000                       for drawing the contour plot.
%      options.PlotXlim =       - A vector specifying the xlim property of
%      'auto'                     the contour plot.
%      options.PlotYlim =       - A vector specifying the ylim property of
%      'auto'                     the contour plot.
%      options.TolBal = 1e-6    - To balance eig(A,B), it is attempted to
%                                 scale |A|.^2+|B|.^2 into a doubly
%                                 stochastic matrix. This happens
%                                 iteratively until the relative
%                                 improvement falls below options.TolBal.
%                                 If solutions are missing, decreasing
%                                 options.TolBal may help.
%      options.TolStep = 1e-2   - The maximum relative step size in the
%                                 Newton-Raphson refinement. Solutions that
%                                 are updated with a larger step size are
%                                 discarded.
%      options.TolSV = 1e-6     - The tolerance for inverting components of
%                                 the Jacobian in the Newton-Raphson
%                                 refinement.
%      options.Univariate =     - True if p and q are to be interpreted as
%      ~isreal(p) || ~isreal(q)   polyanalytic univariate polynomials.
%                                 False if they are bivariate polynomials.
%
%   See also polymin, polymin2, polyval2.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Numerical solution of
%       bivariate and polyanalytic polynomial systems", ESAT-SISTA Internal
%       Report 13-84, KU Leuven, 2013.
%   [2] L. Sorber, M. Van Barel, L. De Lathauwer, "Unconstrained
%       optimization of real functions in complex variables," SIAM J. Opt.,
%       Vol. 22, No. 3, 2012, pp. 879-898.

% Check the options structure.
if nargin < 3, options = struct; end
if ~isfield(options,'MaxIter'),    options.MaxIter  = 8;      end
if ~isfield(options,'Plot'),       options.Plot     = false;  end
if ~isfield(options,'PlotPoints'), options.PlotPoints = 1000; end
if ~isfield(options,'PlotXlim'),   options.PlotXlim = 'auto'; end
if ~isfield(options,'PlotYlim'),   options.PlotYlim = 'auto'; end
if ~isfield(options,'TolBal'),     options.TolBal   = 1e-8;   end
if ~isfield(options,'TolStep'),    options.TolStep  = 1e-2;   end
if ~isfield(options,'TolSV'),      options.TolSV    = 1e-6;   end
if ~isfield(options,'Univariate')
    options.Univariate = ~isreal(p) || ~isreal(q);
end

% Center the exponents around 0 without introducing round-off error.
p = trim(p*2^(-round(median(log2(abs(p(p ~= 0)))))));
q = trim(q*2^(-round(median(log2(abs(q(q ~= 0)))))));

% If bivariate, balance (x,y) and then convert to (z,conj(z)).
if ~options.Univariate
    [bx,by] = balxy(p,q);
    up = bsxfun(@times,bsxfun(@times,p, ...
         bx.^(0:size(p,2)-1)),by.^(0:size(p,1)-1).');
    uq = bsxfun(@times,bsxfun(@times,q, ...
         bx.^(0:size(q,2)-1)),by.^(0:size(q,1)-1).');
	up = bi2uni(up);
    uq = bi2uni(uq);
else
    up = p;
    uq = q;
end

% Balance (z,conj(z)).
dp = size(up)-1;
dq = size(uq)-1;
[bzconjz,bz] = balzz(up,uq);
up = bsxfun(@times,bsxfun(@times,up, ...
     (bzconjz*bz).^(0:dp(2))),bzconjz.^(0:dp(1)).');
uq = bsxfun(@times,bsxfun(@times,uq, ...
     (bzconjz*bz).^(0:dq(2))),bzconjz.^(0:dq(1)).');

% Center the exponents around 0 without introducing round-off error.
up = up*2^(-round(median(log2(abs(up(up ~= 0))))));
uq = uq*2^(-round(median(log2(abs(uq(uq ~= 0))))));

% Generate the Sylvester matrix polynomial in the format [Sd;...;S1;S0].
n = dp(1)+dq(1);
S = zeros(n*(max(dp(2),dq(2))+1),n);
for d = 1:dq(1)
    Sd = [zeros(d-1,dp(2)+1);fliplr(up);zeros(dq(1)-d,dp(2)+1)];
    S(end-numel(Sd)+1:end,d) = Sd(:);
end
for d = 1:dp(1)
    Sd = [zeros(d-1,dq(2)+1);fliplr(uq);zeros(dp(1)-d,dq(2)+1)];
    S(end-numel(Sd)+1:end,dq(1)+d) = Sd(:);
end
Sd = S(1:n,:);
S = -S(n+1:end,:);

% Set up, balance and solve the pencil (A,B) as (Dl*A*Dr,Dl*B*Dr). MATLAB's
% eig(A,B) calls the LAPACK driver routine ZGGEV, which unfortunately does
% not ask ZGGBAL to balance by default.
[Dl,Dr] = balab(S,Sd,n,options);
A = diag(ones(1,size(S,1)-n),n);
A(:,1:n) = S;
B = eye(size(A,1));
B(1:n,1:n) = Sd;
A = bsxfun(@times,bsxfun(@times,A,Dl(:)),Dr(:).');
B = bsxfun(@times,bsxfun(@times,B,Dl(:)),Dr(:).');
z = eig(A,B);

% Undo balancing.
z = z(isfinite(z));
if options.Univariate
    xy = (bzconjz*bz)*z;
else
    xy = [(bzconjz*bz*bx)*real(z) (bzconjz*bz*by)*imag(z)];
end

% Iteratively refine and filter with complex Newton-Raphson.
switch nargout
   case {0,1},xy = newton(p,q,xy,options);
   case 2,   [xy,res] = newton(p,q,xy,options);
   case 3,   [xy,res,dxy] = newton(p,q,xy,options);
end
if nargout > 1
    res = res./[polyval2(abs(p),abs(xy)) polyval2(abs(q),abs(xy))];
end

% Optionally plot contour lines.
if options.Plot && size(xy,1) > 0
    
    % Compute function value.
    if options.Univariate
        x = real(xy); y = imag(xy);
    else
        x = xy(:,1); y = xy(:,2);
    end
    xl = [min(x) max(x)]; yl = [min(y) max(y)];
    xe = 0.15*max(1,xl(2)-xl(1)); ye = 0.15*max(1,yl(2)-yl(1));
    if ischar(options.PlotXlim), xl = [xl(1)-xe xl(2)+xe];
    else xl = options.PlotXlim; end
    if ischar(options.PlotYlim), yl = [yl(1)-ye yl(2)+ye];
    else yl = options.PlotYlim; end
    [X,Y] = meshgrid(linspace(xl(1),xl(2),options.PlotPoints), ...
                     linspace(yl(1),yl(2),options.PlotPoints));
    if options.Univariate
        Zp = reshape(polyval2(p,complex(X(:),Y(:))),size(X));
        Zq = reshape(polyval2(q,complex(X(:),Y(:))),size(X));
    else
        Zp = reshape(polyval2(p,[X(:) Y(:)]),size(X));
        Zq = reshape(polyval2(q,[X(:) Y(:)]),size(X));
    end
    
    % Plot contour lines.
    if ~options.Univariate
        type = 0; % p(x,y), q(x,y)
    elseif all(fliplr(size(p)) == size(q)) && all(all(p == q'))
        type = 1; % p(z,z') = conj(q(z,z'))
    elseif size(p,1) == size(p,2) && size(q,1) == size(q,2) && ...
            ~any(any(p-p')) && ~any(any(q-q'));
        type = 2; % imag(p(z,z')) = imag(q(z,z')) = 0
    else
        type = 3; % p(z,z'), q(z,z')
    end
    if type == 1
        contour(X,Y,real(Zp),[0 0],'b'); hold on;
        contour(X,Y,imag(Zp),[0 0],'r');
    else
        contour(X,Y,real(Zp),[0 0],'b'); hold on;
        contour(X,Y,real(Zq),[0 0],'r');
    end
    if type == 3
        contour(X,Y,imag(Zp),[0 0],'g');
        contour(X,Y,imag(Zq),[0 0],'m');
    end
    
    % Plot solutions.
    plot(x,y,'ok','MarkerSize',4); hold off;
    if options.Univariate
        xlabel('real(z)'); ylabel('imag(z)');
    else
        xlabel('x'); ylabel('y');
    end
    
    % Add legend.
    switch type
        case 0, legend('p(x,y) = 0','q(x,y) = 0','(x*,y*)');
        case 1, legend('real(p(z)) = 0','imag(p(z)) = 0','z*');
        case 2, legend('p(z) = 0','q(z) = 0','z*');
        case 3, legend('real(p(z)) = 0','real(q(z)) = 0', ...
                       'imag(p(z)) = 0','imag(q(z)) = 0','z*');
    end
    
end

end

function [bx,by] = balxy(p,q)
%BALXY Balance p(x,y) and q(x,y) by minimizing the variance in the
% exponents of their anti-diagonals with the transformation x = bx*x and
% y = by*y for some scalars bx and by.

% Set up least-squares system.
p = log2(abs(p)); q = log2(abs(q)); p = p(:); q = q(:);
dxp = repmat(0:size(p,2)-1,size(p,1),1); dxp = dxp(:);
dyp = repmat((0:size(p,1)-1).',1,size(p,2)); dyp = dyp(:);
dxq = repmat(0:size(q,2)-1,size(q,1),1); dxq = dxq(:);
dyq = repmat((0:size(q,1)-1).',1,size(q,2)); dyq = dyq(:);
ifp = isfinite(p); ifq = isfinite(q);
Ab = [dxp(ifp) dyp(ifp) p(ifp); ...
      dxq(ifq) dyq(ifq) q(ifq)];
Ab = bsxfun(@minus,Ab,mean(Ab));
ATA = Ab(:,1:2)'*Ab(:,1:2);
ATb = -Ab(:,1:2)'*Ab(:,3);

% Solve least-squares system for the balancing exponents bxy.
bxy = pinv(ATA,1e-3*norm(ATA,'fro'))*(ATb);
bxy = [floor(bxy) ceil(bxy) [floor(bxy(1));ceil(bxy(2))] ...
       [ceil(bxy(1));floor(bxy(2))]];
res = ATA*bxy-ATb(:,ones(1,size(bxy,2)));
[~,idx] = min(dot(res,res));
bx = 2^bxy(1,idx);
by = 2^bxy(2,idx);

end

function [bzconjz,bz] = balzz(p,q)
%BALZZ Balance p(z,conj(z)) and q(z,conj(z)) by minimizing the variance in
% their exponents with the transformation z = bzconjz*bz*z and conj(z) =
% bzconjz*conj(z) for some scalars bzconjz and bz.
 
% Set up least-squares system.
p = log2(abs(p)); q = log2(abs(q)); p = p(:); q = q(:);
dzp = repmat(0:size(p,2)-1,size(p,1),1); dzp = dzp(:);
dzzp = dzp+repmat((0:size(p,1)-1).',1,size(p,2)); dzzp = dzzp(:);
dzq = repmat(0:size(q,2)-1,size(q,1),1); dzq = dzq(:);
dzzq = dzq+repmat((0:size(q,1)-1).',1,size(q,2)); dzzq = dzzq(:);
ifp = isfinite(p); ifq = isfinite(q);
Ab = [dzzp(ifp) dzp(ifp) p(ifp); ...
      dzzq(ifq) dzq(ifq) q(ifq)];
Ab = bsxfun(@minus,Ab,mean(Ab));
ATA = Ab(:,1:2)'*Ab(:,1:2);
ATb = -Ab(:,1:2)'*Ab(:,3);
 
% Solve least-squares system for the balancing exponents.
bzz = pinv(ATA,1e-3*norm(ATA,'fro'))*(ATb);
bzz = [floor(bzz) ceil(bzz) [floor(bzz(1));ceil(bzz(2))] ...
       [ceil(bzz(1));floor(bzz(2))]];
res = ATA*bzz-ATb(:,ones(1,size(bzz,2)));
[~,idx] = min(dot(res,res));
bzconjz = 2^bzz(1,idx);
bz      = 2^bzz(2,idx);
 
end

function [Dl,Dr] = balab(S,Sd,n,options)
%BALAB Balance eig(A,B) as eig(Dl*A*Dr,Dl*B*Dr) by scaling |A|.^2+|B|.^2
% into a doubly stochastic matrix. I.e., its row and column sums should
% equal 1. The scaling matrices Dl and Dr are rounded to the nearest powers
% of 2.

% For the pencil (A,B), compute (part of) the matrix M = |A|.^2+|B|.^2.
M = S.*conj(S);
M(1:size(Sd,1),1:size(Sd,2)) = M(1:size(Sd,1),1:size(Sd,2))+Sd.*conj(Sd);

% Alternate between updates for Dl and Dr.
xk = ones(2*size(M,1),1);
Axk = Ax(xk,false);
dy = max(abs(xk(1:end/2).*Axk-1));
for i = 1:2*size(xk,1)
    xk(1:end/2) = 1./Axk;
    Axk = Ax(xk,true);
    xk(end/2+1:end) = 1./Axk;
    Axk = Ax(xk,false);
    dy1 = dy;
    dy = max(abs(xk(1:end/2).*Axk-1));
    if 1-dy/dy1 <= options.TolBal, break; end
end

% Round Dl and Dr.
Dl = 2.^round(0.5*log2(abs(xk(1:end/2))));
Dr = 2.^round(0.5*log2(abs(xk(end/2+1:end))));

function y = Ax(x,transp)
    if ~transp
        % Computes y = [0 M]*x.
        x2 = x(end/2+1:end);
        y = M*x2(1:size(S,2));
        y(1:end-n) = y(1:end-n)+x2(n+1:end);
        y(n+1:end) = y(n+1:end)+x2(n+1:end);
    else
        % Computes y = [M.' 0]*x.
        x1 = x(1:end/2);
        y = zeros(size(S,1),1);
        y(1:n) = x1.'*M;
        y(n+1:end) = y(n+1:end)+x1(1:end-n)+x1(n+1:end);
    end
end

end

function up = bi2uni(p)
%BI2UNI Convert a polynomial p in two real variables (x,y) to a real
% polyanalytic polynomial in one complex variable (z,conj(z)), where
% z := x+1i*y, using a variant of Pascal's triangle.
	
% Convert polynomial assuming p(end) ~= 0, i.e. maximal total degree.
n = sum(size(p))-1;
up = zeros(n,n);
C = zeros(n,n); C(1,1) = 1;
for j = 1:size(p,2)
    if j > 1
        % Compute Pascal's triangle up to the j-th antidiagonal.
        deg = j-1;
        iprv = deg:n-1:n*(deg-1)+1;
        inxt = deg+1:n-1:n*deg+1;
        C(inxt) = [0 C(iprv)]+[C(iprv) 0];
    end
    for i = 1:size(p,1)
        % Compute temporary Pascal's triangle variant starting from the
        % (j+1)-th diagonal (using (z-conj(z)) instead of (z+conj(z))).
        deg = i+j-2;
        inxt = deg+1:n-1:n*deg+1;
        if i > 1
            iprv = deg:n-1:n*(deg-1)+1;
            C(inxt) = [0 C(iprv)]-[C(iprv) 0];
        end
        if p(i,j) ~= 0
            up(inxt) = up(inxt)+p(i,j)*(0.5^(i+j-2)*(-1i)^(i-1))*C(inxt);
        end
    end
end

% Trim the polynomial, this is necessary for cases where the total
% degree of p is less than maximal.
up = trim(up);

end

function p = trim(p)
%TRIM Removes trailing zero columns and rows from a polynomial.

for i = size(p,1):-1:1, if any(p(i,:) ~= 0), break; end; end
for j = size(p,2):-1:1, if any(p(:,j) ~= 0), break; end; end
p = p(1:i,1:j);

end

function [z,F,J] = newton(p,q,z0,options)
%NEWTON Locate roots of the polynomial system p = q = 0 with a complex
% Newton-Raphson algorithm. The algorithm terminates if (a) it can not take
% a step that improves the residual, (b) the step size is larger than the
% specified tolerance or (c) the maximum number of iterations has been
% reached. Solutions for which the step size is larger than the specified
% tolerance are removed entirely.

% Newton-Raphson.
dpdz = trim(p(:,2:end)*diag(1:size(p,2)-1));
dqdz = trim(q(:,2:end)*diag(1:size(q,2)-1));
dpdconjz = diag(1:size(p,1)-1)*p(2:end,:);
dqdconjz = diag(1:size(q,1)-1)*q(2:end,:);
F = @(z)[polyval2(p,z) polyval2(q,z)];
J = @(z)[polyval2(dpdz,z) polyval2(dqdz,z) ...
         polyval2(dpdconjz,z) polyval2(dqdconjz,z)];
K = [1 1i; 1 -1i];
z = z0;
stop = false(size(z,1),1);
for i = 1:options.MaxIter

    % Compute Jacobian.
    Fz = F(z(~stop,:));
    Jz = J(z(~stop,:));
    
    % Compute the Newton-Raphson step.
    dz = zeros(sum(~stop),size(z,2));
    if options.Univariate
        for j = 1:sum(~stop)
            A = reshape(Jz(j,:),2,2)*K;
            A = [real(A(1,:));imag(A(1,:));real(A(2,:));imag(A(2,:))];
            b = -[real(Fz(j,1));imag(Fz(j,1));real(Fz(j,2));imag(Fz(j,2))];
            dxy = [0 0];
            if ~all(isfinite(A(:))), continue; end
            [u,s,v] = svd(A);
            sv = [s(1) s(end)/s(1)] > options.TolSV;
            if sv(1) && sv(2)
                dxy = v*diag(1./diag(s))*(u'*b);
            elseif sv(1)
                dxy = v(:,1)/s(1)*(u(:,1)'*b);
            end
            dz(j) = dxy(1)+dxy(2)*1i;
        end
        smallstep = all((abs([real(z(~stop)) imag(z(~stop))]) > 1 & ...
            abs([real(dz) imag(dz)]) <= ...
            options.TolStep*abs([real(z(~stop)) imag(z(~stop))])) | ...
            (abs([real(z(~stop)) imag(z(~stop))]) <= 1 & ...
            abs([real(dz) imag(dz)]) <= options.TolStep),2);
    else
        for j = 1:sum(~stop)
            A = reshape(Jz(j,:),2,2);
            b = -Fz(j,:).';
            if ~all(isfinite(A(:))), continue; end
            [u,s,v] = svd(A);
            sv = [s(1) s(end)/s(1)] > options.TolSV;
            if sv(1) && sv(2)
                dz(j,:) = v*diag(1./diag(s))*(u'*b);
            elseif sv(1)
                dz(j,:) = v(:,1)/s(1)*(u(:,1)'*b);
            end
        end
        smallstep = all((abs(z(~stop,:)) > 1 & abs(dz) <= ...
            options.TolStep*abs(z(~stop,:))) | ...
            (abs(z(~stop,:)) <= 1 & abs(dz) <= options.TolStep),2);
    end
    
    % Discard unlikely candidates.
    if ~all(smallstep)
        dz = dz(smallstep,:);
        keep = stop;
        keep(~stop) = smallstep;
        stop = stop(keep);
        z = z(keep,:);
    end
    if all(stop), break; end
    z(~stop,:) = z(~stop,:)+dz;
    
end

% Compute residuals and Jacobian.
if nargout > 1, F = F(z); end
if nargout > 2, J = reshape(J(z).',2,2,[]); end

end
