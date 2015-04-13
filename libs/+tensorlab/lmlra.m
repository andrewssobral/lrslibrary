function [U,S,output] = lmlra(T,size_core,S0,options)
%LMLRA Low multilinear rank approximation.
%   [U,S,output] = lmlra(T,size_core) computes the factor matrices U{1},
%   ..., U{N} and core tensor S of dimensions size_core belonging to a low
%   multilinear rank approximation of the N-th order tensor T.
%
%   [U,S,output] = lmlra(T,U0,S0) allows the user to provide initial factor
%   matrices U0 and core tensor S0, to be used by the main algorithm.
%   Factor matrices equal to the empty matrix [] will be initialized using
%   the chosen initialization method.
%
%   The structure output contains the output of the selected algorithms:
%
%      output.Initialization - The output of the initialization step.
%      output.Algorithm      - The output of the main algorithm.
%
%   lmlra(T,size_core,options) and lmlra(T,U0,S0,options) allow the user to
%   choose which algorithms will be used in the different steps and also
%   set parameters for those algorithms:
%
%      Use options.Initialization = [{@lmlra_aca}|@lmlra_rnd|@mlsvd] to
%      choose the initialization method. By default, an adaptive cross
%      approximation algorithm is used. The variable
%      options.InitializationOptions will be passed to the chosen 
%      initialization method as third argument.
%
%      Use options.Algorithm = [@lmlra_hooi|@lmlra_minf|{@lmlra_nls}| ...
%      lmlra3_dgn|lmlra3_rtr] to choose the main algorithm. The structure
%      options.AlgorithmOptions will be passed to the chosen algorithm.
%
%   Further options are:
%
%      options.Display = false - Set to true to enable printing output
%                                information to the command line.
%
%   See also mlrankest, lmlragen, lmlraerr.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Gather some input information.
T = fmt(T);
N = length(size_core);
if nargin < 3, S0 = []; end
if iscell(size_core) || (nargin >= 3 && isnumeric(S0) && ~isempty(S0))
    U = size_core;
    S = S0;
    size_core = cellfun('size',U,2);
end

% Check the options structure.
isfunc = @(f)isa(f,'function_handle');
xsfunc = @(f)isfunc(f)&&exist(func2str(f),'file');
if nargin < 4
    if nargin == 3 && isstruct(S0), options = S0;
    else options = struct; end
end
if ~isfield(options,'Initialization')
    funcs = {@lmlra_aca,@lmlra_rnd,@mlsvd};
    options.Initialization = funcs{find(cellfun(xsfunc,funcs),1)};
end
if length(size_core) == 1, options.Initialization = @mlsvd; end
if ~isfield(options,'Algorithm')
    funcs = {@lmlra_nls,@lmlra_minf,@lmlra_hooi,@lmlra3_rtr,@lmlra3_dgn};
    options.Algorithm = funcs{find(cellfun(xsfunc,funcs),1)};
end
if ~isfield(options,'Display'), options.Display = false; end
if ~options.Display, print = @(varargin)true; else print = @fprintf; end

% Step 1: initialize the factor matrices unless they were all provided by
% the user.
print('Step 1: Initialization ');
if ~exist('U','var'), U = cell(1,N); end
Uempty = cellfun(@isempty,U);
if any(Uempty) || isempty(S0)
    if isfunc(options.Initialization)
        print('is %s... ',func2str(options.Initialization));
    end
    if ~xsfunc(options.Initialization)
        error('lmlra:Initialization','Not a valid initialization.');
    end
else
    options.Initialization = false;
    print('is manual... ');
end
if isfunc(options.Initialization)
    % Generate initial factor matrices.
    if ~isfield(options,'InitializationOptions')
        options.InitializationOptions = struct;
    end
    [U0,S0,sv] = options.Initialization(T,size_core,...
                 options.InitializationOptions);
    if isstruct(sv)
        output.Initialization = sv;
    else
        output.Initialization.sv = sv;
    end
    % Fill empty factor matrices.
    for n = 1:sum(Uempty), U{n} = U0{n}; end
    if ~isempty(S0), S = S0; end
    output.Initialization.Name = func2str(options.Initialization);
else
    output.Initialization.Name = 'manual';
end
output.Initialization.relerr = frob(lmlrares(T,U,S))/frob(T);
print('relative error = %.6g.\n',output.Initialization.relerr);

% Step 2: run the main LMLRA algorithm.
if xsfunc(options.Algorithm)
    print('Step 2: Algorithm is %s... ',func2str(options.Algorithm));
    if ~isfield(options,'AlgorithmOptions')
        options.AlgorithmOptions = struct;
    end
    [U,S,output.Algorithm] = options.Algorithm(T,U,S,...
                             options.AlgorithmOptions);
    output.Algorithm.Name = func2str(options.Algorithm);
else
    error('lmlra:Algorithm','Not a valid algorithm.');
end
if isfield(output.Algorithm,'iterations')
    print('iterations = %i, ',output.Algorithm.iterations);
end
output.Algorithm.relerr = frob(lmlrares(T,U,S))/frob(T);
print('relative error = %.6g.\n',output.Algorithm.relerr);

% Format output.
fn = fieldnames(output);
for f = 1:length(fn)
    output.(fn{f}) = orderfields(output.(fn{f}));
end
