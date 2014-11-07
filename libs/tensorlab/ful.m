function T = ful(T)
%FUL Convert formatted data set to an array.
%   T = ful(T) converts the formatted data set T into a MATLAB array. If T
%   is incomplete, unknown entries are represented by the value NaN.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.

if isnumeric(T), return; end
if ~isstruct(T)
    error('ful:T','T must be a data set formatted by fmt.');
end

val = T.val;
size_tens = T.size;
if 8*prod(size_tens) > 8e9
    error('ful:size','T is too large to be stored as an array.');
end
if ~isfield(T,'ind')
    ind = sub2ind(size_tens,T.sub{:});
else
    ind = T.ind;
end

if isfield(T,'sparse') && T.sparse
    T = zeros(size_tens);
elseif isfield(T,'incomplete') && T.incomplete
    T = nan(size_tens);
else
    error('ful')
end
T(ind) = val;
