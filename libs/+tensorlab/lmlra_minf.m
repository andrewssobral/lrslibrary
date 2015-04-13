function [U,S,output] = lmlra_minf(T,U0,S0,options)
%LMLRA_MINF LMLRA by unconstrained nonlinear optimization.
%   [U,S,output] = lmlra_minf(T,U0,S0) computes the factor matrices U{1},
%   ..., U{N} and core tensor S belonging to a low multilinear rank
%   approximation of the N-th order tensor T by minimizing 
%   0.5*frob(T-lmlragen(U,S))^2. Each term U{r} is a cell array of N factor
%   matrices U{r}{n}, followed by a core tensor U{r}{N+1}. The algorithm is
%   initialized with the factor matrices U0{n} and core tensor S0. The
%   structure output returns additional information:
%
%      output.Name  - The name of the selected algorithm.
%      output.<...> - The output of the selected algorithm.
%
%   lmlra_minf(T,U0,S0,options) may be used to set the following options:
%
%      options.Algorithm =     - The desired optimization method.
%      [{@minf_lbfgsdl}|...
%       @minf_lbfgs|@minf_ncg]
%      options.<...>           - Parameters passed to the selected method,
%                                e.g., options.TolFun, options.TolX and
%                                options.LineSearchOptions. See also help
%                                [options.Algorithm].
%
%   See also lmlra_nls.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.
%   [2] L. Sorber, M. Van Barel, L. De Lathauwer, "Unconstrained
%       optimization of real functions in complex variables," SIAM J. Opt.,
%       Vol. 22, No. 3, 2012, pp. 879-898.

% Forward problem to btd_minf as a one-term BTD.
if nargin < 4, options = struct; end
[U,output] = btd_minf(T,{[U0(:).',S0]},options);
S = U{1}{end};
U = U{1}(1:end-1);
