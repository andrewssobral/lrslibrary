function [U,output] = btd_rnd(size_tens,size_core,options)
%BTD_RND Pseudorandom initialization for BTD.
%   U = btd_rnd(size_tens,size_core) generates R pseudorandom terms U{r} to
%   initialize algorithms that compute a block term decomposition of an
%   N-th order tensor. Each term U{r} is a cell array of N unitary factor
%   matrices U{r}{n}, followed by a core tensor U{r}{N+1} of size
%   size_core{r}.
%
%   btd_rnd(T,size_core) is shorthand for btd_rnd(size(T),size_core) if
%   T is real. If T is complex, then by default the terms will be generated
%   using pseudorandom complex numbers as well (cf. options).
%
%   btd_rnd(size_tens,size_core,options) and btd_rnd(T,size_core, ...
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
%   See also cpd_rnd, lmlra_rnd, btdgen.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

if nargin < 3, options = struct; end
if ~iscell(size_core), size_core = {size_core}; end
function U = capture(size_core)
    [U,S] = lmlra_rnd(size_tens,size_core,options);
    U{end+1} = S;
end
U = cellfun(@capture,size_core,'UniformOutput',false);
output = struct;

end
