function [U,output] = cpd_rnd(size_tens,R,options)
%CPD_RND Pseudorandom initialization for CPD.
%   U = cpd_rnd(size_tens,R) generates pseudorandom factor matrices U{1},
%   ..., U{N} of dimensions size_tens(n)-by-R that can be used to
%   initialize algorithms that compute the CPD of an N-th order tensor.
%
%   cpd_rnd(T,R) is shorthand for cpd_rnd(size(T),R) if T is real. If T is
%   complex, then U{n} will be generated as complex matrices by default
%   (cf. options).
%
%   cpd_rnd(size_tens,R,options) and cpd_rnd(T,R,options) may be used to
%   set the following options:
%
%      options.Real =        - The type of random number generator used to
%      [{@randn}|@rand|0]      generate the real part of each factor
%                              matrix. If 0, there is no real part.
%      options.Imag =        - The type of random number generator used to
%      [@randn|@rand|0|...     generate the imaginary part of each factor
%       {'auto'}]              matrix. If 0, there is no imaginary part.
%                              On 'auto', options.Imag is 0 unless the
%                              first argument is a complex tensor T, in
%                              which case it is equal to options.Real.
%      options.Orth =        - If true, the generated factor matrices are
%      [true|false|{'auto'}]   orthogonalized using a QR factorization.
%                              On 'auto', options.Orth is false if the
%                              first argument is a vector size_tens, and
%                              true if the first argument is a tensor T.
%
%   See also btd_rnd, lmlra_rnd, cpdgen.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Process input.
isSizeVector = isnumeric(size_tens) && isvector(size_tens);
if ~isSizeVector
    T = size_tens;
    if isstruct(T), size_tens = T.size;
    else size_tens = size(T); end
end
N = length(size_tens);

% Check the options structure.
isfunc = @(f)isa(f,'function_handle');
if nargin < 3, options = struct; end
if ~isfield(options,'Real'), options.Real = @randn; end
if ~isfunc(options.Real), options.Real = @zeros; end
if ~isfield(options,'Imag'), options.Imag = 'auto'; end
if ischar(options.Imag) && strcmpi(options.Imag,'auto')
    if ~isSizeVector && ((isstruct(T) && ~isreal(T.val)) || ...
            (isnumeric(T) && ~isreal(T)))
        options.Imag = options.Real;
    else
        options.Imag = 0;
    end
end
if ~isfield(options,'Orth'), options.Orth = 'auto'; end
if ischar(options.Orth) && strcmpi(options.Orth,'auto')
	options.Orth = ~isSizeVector;
end

% Generate factor matrices.
U = arrayfun(@(n)options.Real(size_tens(n),R),1:N,'UniformOutput',0);
if isfunc(options.Imag)
    Ui = arrayfun(@(n)options.Imag(size_tens(n),R),1:N,'UniformOutput',0);
    U = cellfun(@(ur,ui)ur+ui*1i,U,Ui,'UniformOutput',0);
end
for n = 1:N*options.Orth
    if size(U{n},1) >= size(U{n},2), [U{n},~] = qr(U{n},0);
    else [Q,~] = qr(U{n}.',0); U{n} = Q.'; end
end
output = struct;
