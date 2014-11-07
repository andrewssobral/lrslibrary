function [alpha,output] = cpd_els(T,U,dU,state,options)
%CPD_ELS CPD exact line search.
%   [alpha,output] = cpd_els(T,U,dU) computes the real or complex line
%   search parameter alpha that minimizes the objective function f(alpha)
%   := 0.5*frob(T-cpdgen(U+alpha*dU))^2, in which T is an N-th order
%   tensor, U is a cell array containing N factor matrices U{n}, and dU is
%   a cell array containing N step matrices dU{n}.
%
%   If options.Scale is true, the objective function f(alpha) := 0.5* ...
%   frob(T-cpdgen(alpha(2)*(U+alpha(1)*dU)))^2 is minimized for both
%   parameters alpha(1) and alpha(2).
%
%   Note that for incomplete tensors, the objective function is partially
%   approximated and the alpha is no longer a global minimizer in general.
%
%   The structure output returns additional information:
%
%      output.fval    - The value of the objective function f at alpha.
%      output.relgain - The improvement in objective function w.r.t.
%                       alpha = 0, relative to the improvement of the step
%                       alpha = 1 w.r.t. alpha = 0.
%
%   cpd_els(T,U,dU,[],options) may be used to set the following options:
%
%      options.IsReal =      - If true, alpha is restricted to the real
%      isreal(T)               domain. Otherwise, alpha may be complex
%                              (e.g., in case of a complex CPD).
%      options.Scale = false - If true, computes an optimal scaling
%                              parameter alpha(2) in conjunction with an
%                              optimal line search parameter alpha(1).
%      options.SolverOptions - A struct passed to the the polynomial or
%                              rational minimizer. See polymin, polymin2, 
%                              ratmin and ratmin2 for details.
%
%   See also cpd_aels, cpd_eps.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, I. Domanov, M. Van Barel, L. De Lathauwer, "Exact Line
%       and Plane Search for Tensor Optimization", ESAT-SISTA Internal
%       Report 13-02, KU Leuven, 2013.
%   [2] M. Rajih, P. Comon, R. Harshman, "Enhanced line search: A novel
%       method to accelerate PARAFAC," SIAM J. Matrix Anal. Appl., Vol. 30,
%       2008, pp. 1148-1171.

% Check the factor matrices U and dU.
N = length(U);
R = size(U{1},2);
if length(U) ~= N || length(dU) ~= N
    error('cpd_eps:U','length((d)U) should equal ndims(T).');
end
if any(cellfun('size',U,2) ~= R) || any(cellfun('size',dU,2) ~= R)
    error('cpd_els:U','size((d)U{n},2) should be the same for all n.');
end

% Check the options structure.
if nargin < 5, options = struct; end
if ~isfield(options,'IsReal')
    options.IsReal = (isstruct(T) && isreal(T.val)) || ....
        (~isstruct(T) && isreal(T));
end
if ~isfield(options,'Scale'), options.Scale = false; end
if ~isfield(options,'SolverOptions'), options.SolverOptions = struct; end
if ~isfield(options.SolverOptions,'TolBal')
    options.SolverOptions.TolBal = 1e-2;
end

% Check the state structure.
if nargin < 4 || ~isstruct(state), state = struct; end
if ~isfield(state,'fval') || ...
   (options.Scale && (~isfield(state,'K') || ~isfield(state,'K1')))
    if ~isfield(state,'fval')
        D = cpdres(T,U);
        state.fval = 0.5*(D(:)'*D(:));
    end
    if options.Scale && (~isfield(state,'K') || ~isfield(state,'K1'))
        state.K1 = 1;
        state.K{state.K1} = mtkrprod(T,U,1);
    end
end
if ~isfield(state,'T2')
    state.T2 = frob(T,'squared');
end

% Generate indices of combinations of terms.
Uall = [U(:).'; dU(:).'];
I = cell2mat(arrayfun(@(i)bitget(i,1:N),(0:2^N-1)','UniformOutput',false));
[deg,idx] = sort(sum(I,2));
I = I(idx,:)+1;

% Compute the polynomial p = <T,cpdgen(U+z*dU)> if options.Scale is true.
% Otherwise compute p = <cpdgen(U)-T,cpdgen(U+z*dU)> so that the objective
% function f is more accurate later on. This costs about 2^N-1 function
% evaluations (but can be done less accurately in N evaluations by
% interpolating on the unit circle).
p = zeros(1,N+1);
if options.Scale
    p(1) = sum(dot(state.K{state.K1},U{state.K1}));
else
    T = cpdres(T,U);
end
for d = 1:N
    for i = find(deg.' == d)
        tmp = conj(sum(mtkrprod(T,Uall(I(i,1:end)+2*(0:N-1)),0)));
        p(d+1) = p(d+1)+tmp;
    end
end

% Build a cell array of all products (d)U{n}'*(d)U{n}.
UHU = cell(1,N);
for n = 1:N
    tmp = cat(2,U{n},dU{n});
    UHU{n} = tmp'*tmp;
end

% Compute the polynomial q = ||cpdgen(U+z*dU)||^2.
q = zeros(N+1,N+1);
for j = 1+~options.Scale:N+1
    for i = 1+~options.Scale:j
        for n = 1:N
            col = bsxfun(@plus,1:R,R*(I(deg == j-1,n)-1)).';
            row = bsxfun(@plus,1:R,R*(I(deg == i-1,n)-1)).';
            UHUn = UHU{n};
            if n == 1, W = UHUn(row,col); else W = W.*UHUn(row,col); end
        end
        q(i,j) = sum(W(:));
    end
end
q = triu(q,1)+diag(real(diag(q)))+triu(q,1)';

% Compute objective function f = 0.5*[s*s'*q - s*p - s'*conj(p) + ||T||^2].
f = q;
if options.Scale
    f(1,:) = f(1,:)-p;
    f(:,1) = f(1,:)';
else
    f(1,:) = p;
	f(:,1) = p';
end
f = 0.5*f;
f(1) = state.fval(end);

% Solve (scaled) line search.
if options.Scale
    
    if options.IsReal
        
        % Compute real alpha = argmin -real(p)^2/q.
        pr = fliplr(real(p));
        qr = zeros(1,N+1);
        for d = -N:N, qr(d+N+1) = sum(real(diag(fliplr(q),d))); end
        [z,v] = ratmin(-conv(pr,pr),qr,options.SolverOptions);
        v(state.T2+v < -1e-6*state.T2) = inf;
        [output.fval,idx] = min(v);
        output.fval = 0.5*(state.T2+output.fval);
        alpha = z(idx);
        warn = isempty(alpha);
        if ~warn
            alpha(2) = nthroot(polyval(pr,alpha)/polyval(qr,alpha),N);
        end
        
    else
        
        % Compute complex alpha = argmin -p*conj(p)/q.
        options.SolverOptions.Univariate = true;
        [z,v] = ratmin2(-p'*p,q,options.SolverOptions);
        v(state.T2+v < -1e-6*state.T2) = inf;
        [output.fval,idx] = min(v);
        output.fval = 0.5*(state.T2+output.fval);
        alpha = z(idx);
        warn = isempty(alpha);
        if ~warn
            alpha(2) = (polyval(fliplr(p),alpha)'/polyval2(q,alpha))^(1/N);
        end
        
    end
    
else
    
    if options.IsReal
        
        % Compute real alpha = argmin f.
        fr = zeros(1,N+1);
        for d = -N:N, fr(d+N+1) = sum(real(diag(fliplr(f),d))); end
        [z,v] = polymin(fr,options.SolverOptions);
        v(v < -1e-6*state.T2) = inf;
        [output.fval,idx] = min(v);
        alpha = z(idx);
        warn = isempty(alpha);
        
    else
        
        % Compute complex alpha = argmin f.
        options.SolverOptions.Univariate = true;
        [z,v] = polymin2(f,options.SolverOptions);
        v(v < -1e-6*state.T2) = inf;
        [output.fval,idx] = min(v);
        alpha = z(idx);
        warn = isempty(alpha);
        
    end
    
end

% Compute improvement compared to a step of alpha = 1.
if isempty(output.fval) || state.fval(end) < 0 || ...
   log10(state.fval(end)) <= log10(state.T2)-16+2.5
    output = rmfield(output,'fval');
    output.relgain = nan;
else
    fval0 = f(1); fval1 = real(sum(f(:)));
    output.relgain = (fval0-output.fval)/max(0,fval0-fval1);
end
if warn || (isfield(output,'fval') && output.fval > f(1))
    if options.Scale, alpha = [1 1];
    else alpha = 1; end
    output = struct;
end
