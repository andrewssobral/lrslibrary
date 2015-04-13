function [U,S,output] = lmlra_rnd(size_tens,size_core,options)
%LMLRA_RND Pseudorandom initialization for LMLRA.
%   [U,S] = lmlra_rnd(size_tens,size_core) generates pseudorandom
%   unitary factor matrices U{1}, ..., U{N} of dimensions size_tens(n)-by-
%   size_core(n) and a core tensor of dimensions size_core that can be used
%   to initialize algorithms that compute a low multilinear rank
%   approximation of an N-th order tensor.
%
%   lmlra_rnd(T,size_core) is shorthand for lmlra_rnd(size(T),size_core) if
%   T is real. If T is complex, then by default S and U{n} will be
%   generated as a complex tensor and matrices, respectively (cf. options).
%
%   lmlra_rnd(size_tens,size_core,options) and lmlra_rnd(T,size_core, ...
%   options) may be used to set the following options:
%
%      options.Real =         - The type of random number generator used to
%      [{@randn}|@rand|0]       generate the real part of the core tensor S
%                               and matrices U{n}. If 0, there is no real
%                               part.
%      options.Imag =         - The type of random number generator used to
%      [@randn|@rand|0|...      generate the imaginary part of the core
%       {'auto'}]               tensor S and matrices U{n}. If 0, there is
%                               no imaginary part. On 'auto', options.Imag
%                               is 0 unless the first argument is a complex
%                               tensor T, in which case it is equal to
%                               options.Real.
%      options.Unitary = true - Set equal to false to generate the factor
%                               matrices U{n} without converting them to
%                               unitary matrices afterwards.
%
%   See also btd_rnd, cpd_rnd, lmlragen.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Process input.
isSizeVector = any(size(size_tens) == numel(size_tens));
if ~isSizeVector, T = size_tens; size_tens = size(size_tens); end
N = length(size_tens);
if N ~= length(size_core)
    error('lmlra_rnd:size_core', ...
         ['length(size_core) should be equal to ndims(T) or ' ...
          'length(size_tens).']);
end

% Check the options structure.
isfunc = @(f)isa(f,'function_handle');
if nargin < 3, options = struct; end
if ~isfield(options,'Real'), options.Real = @randn; end
if ~isfunc(options.Real), options.Real = @zeros; end
if ~isfield(options,'Imag'), options.Imag = 'auto'; end
if ~isfield(options,'Unitary'), options.Unitary = true; end
if ischar(options.Imag) && strcmpi(options.Imag,'auto')
    if ~isSizeVector && ~isreal(T)
        options.Imag = options.Real;
    else
        options.Imag = 0;
    end
end

% Generate factor matrices and core tensor.
U = arrayfun(@(n)options.Real(size_tens(n),size_core(n)),1:N, ...
        'UniformOutput',0);
if isfunc(options.Imag)
    Ui = arrayfun(@(n)options.Imag(size_tens(n),size_core(n)),1:N, ...
            'UniformOutput',0);
    U = cellfun(@(ur,ui)ur+ui*1i,U,Ui,'UniformOutput',0);
end
if options.Unitary
    for n = 1:N
        if size(U{n},1) >= size(U{n},2), [U{n},~] = qr(U{n},0);
        else [Q,~] = qr(U{n}.',0); U{n} = Q.'; end
    end
end
S = options.Real(size_core(:).');
if isfunc(options.Imag), S = S+1i*options.Imag(size_core(:).'); end
output = struct;
