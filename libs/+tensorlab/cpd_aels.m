function [alpha,output] = cpd_aels(T,U,dU,state,~)
%CPD_AELS CPD approximate enhanced line search.
%   [alpha,output] = cpd_aels(T,U,dU,state) computes the real line search
%   parameter alpha that approximately minimizes the objective function
%   f(alpha) = 0.5*frob(T-cpdgen(U+alpha*dU))^2, in which T is an N-th
%   order tensor, U is a cell array containing N factor matrices U{n}, and
%   dU is a cell array containing N step matrices dU{n}. The method
%   minimizes the polynomial interpolating f at alpha = [0 .5 1], using the
%   equivalent cost of only one objective function evaluation, if possible.
%
%   See also cpd_els.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Check the factor matrices U and dU.
N = length(U);
R = size(U{1},2);
if any(cellfun('size',U,2) ~= R) || any(cellfun('size',dU,2) ~= R)
    error('cpd_aels:U','size((d)U{n},2) should be the same for all n.');
end

% Check the state structure.
if nargin < 4 || ~isstruct(state), state = struct; end
if ~isfield(state,'fval')
    error('cpd_aels:state','State must contain the field "fval".');
end
if ~isfield(state,'fval1')
    fastFval = isfield(state,'K') && isfield(state,'KN') && ...
               isfield(state,'T2') && isfield(state,'UHU');
end

% Compute objective function value at step = [0 .5 1].
fval0 = state.fval(end);
if isfield(state,'fval1'), fval1 = state.fval1;
else
    % The function value at step = 1 is not supplied.
    if fastFval && log10(state.fval(end)) > log10(state.T2)-16+2.5
        % Use a cheap expression for fval1.
        fval1 = abs(0.5*(state.T2+sum(sum(real(prod(state.UHU,3)))))- ...
                    real(sum(dot(state.K{state.KN}, ...
                    U{state.KN}+dU{state.KN}))));
    else
        % Evaluate objective function.
        A = cellfun(@(u,v)u+v,U,dU,'UniformOutput',false);
        D = cpdres(T,A);
        fval1 = 0.5*(D(:)'*D(:));
    end
end
t = 0.5;
A = cellfun(@(u,v)u+t*v,U,dU,'UniformOutput',false);
D = cpdres(T,A);
fvalt = 0.5*(D(:)'*D(:));

% Compute the minimizer of the quadratic interpolating polynomial.
abc = [0 0 1; 1 1 1; t*t t 1]\[fval0;fval1;fvalt];
if abc(1) >= 0
    alpha = -abc(2)/(2*abc(1));
else
    alpha = 1;
end
output = struct;
