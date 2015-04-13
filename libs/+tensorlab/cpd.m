function [U,output] = cpd(T,R,options)
%CPD Canonical polyadic decomposition.
%   [U,output] = cpd(T,R) computes the factor matrices U{1}, ..., U{N}
%   belonging to a canonical polyadic decomposition of the N-th order
%   tensor T in R rank-one terms.
%
%   [U,output] = cpd(T,U0) allows the user to provide initial factor
%   matrices U0, to be used by the main algorithm. By providing an
%   initialization, the compression step will be skipped. Factor matrices
%   equal to the empty matrix [] will be initialized using the chosen
%   initialization method.
%
%   The structure output contains the output of the selected algorithms:
%
%      output.Compression    - The output of the compression step.
%      output.Initialization - The output of the initialization step.
%      output.Algorithm      - The output of the main algorithm.
%      output.Refinement     - The output of the the refinement step after
%                              decompression (if applicable).
%
%   cpd(T,R,options) and cpd(T,U0,options) allow the user to choose which
%   algorithms will be used in the different steps and also set parameters
%   for those algorithms:
%
%      Use options.Compression = [{'auto'}|true|false] to select whether or
%      not the tensor is first compressed to a tensor of dimensions
%      min(size(T),R*ones(1,N)) with adaptive cross approximation. If one
%      or more initial factor matrices are provided compression is skipped.
%      On 'auto', compression is only executed if it is expected to lead to
%      a significant reduction in computational complexity.
%
%      Use options.Initialization = [{'auto'}|@cpd_rnd|@cpd_gevd] to
%      choose the initialization method. The structure
%      options.InitializationOptions will be passed to the chosen
%      initialization method. On 'auto', @cpd_gevd is used when possible,
%      otherwise @cpd_rnd.
%
%      Use options.Algorithm = [@cpd_als|@cpd_minf|{@cpd_nls}|@cpd3_sd|...
%      @cpd3_sgsd|cpdnn_nlsb|cpds_minf|cpds_nls] to choose the main
%      algorithm. The structure options.AlgorithmOptions will be passed to
%      the chosen algorithm.
%
%      Use options.Refinement = [@cpd_als|@cpd_minf|{@cpd_nls}|@cpd3_sd|...
%      @cpd3_sgsd|cpdnn_nlsb|false] to choose the algorithm used to refine
%      the solution after decompression, if applicable. By default,
%      options.Refinement is equal to options.Algorithm. The structure
%      options.RefinementOptions will be passed to the chosen refinement
%      method. Parameters not set in options.RefinementOptions will be
%      copied from options.AlgorithmOptions where possible.
%
%   Further options are:
%
%      options.Display = false - Set to true to enable printing output
%                                information to the command line.
%
%   See also rankest, cpdgen, cpderr.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Check the tensor T.
T = fmt(T);
if isstruct(T), size_tens = T.size;
else size_tens = size(T); end
if iscell(R)
    U = R;
    R = max(cellfun('size',U,2));
    N = length(U);
else
    if isstruct(T), N = length(T.size);
    else N = ndims(T); end
end

% Check the options structure.
isfunc = @(f)isa(f,'function_handle');
xsfunc = @(f)isfunc(f)&&exist(func2str(f),'file');
if nargin < 3, options = struct; end
if ~isfield(options,'Compression')
    options.Compression = 'auto';
end
if ~isfield(options,'Initialization')
    options.Initialization = 'auto';
end
if ~isfield(options,'Algorithm')
    funcs = {@cpd_nls,@cpd_als,@cpd_minf,@cpd3_sd,@cpd3_sgsd};
    options.Algorithm = funcs{find(cellfun(xsfunc,funcs),1)};
end
if ~isfield(options,'Refinement')
    options.Refinement = options.Algorithm;
end
if ~isfield(options,'Display'), options.Display = false; end
if ~options.Display, print = @(varargin)true; else print = @fprintf; end

% STEP 1: compress the tensor if it was requested, and no initial factor
% matrices were provided by the user.
print('Step 1: Compression ');
if ~exist('U','var') && R > 1
    size_core = min(size_tens,R*ones(1,N));
    ratio = prod(size_core)/prod(size_tens);
    if ischar(options.Compression) && strcmpi(options.Compression,'auto')
        if ratio < 0.5
            options.Compression = true;
            print('is lmlra_aca (compression ratio %.6g < 0.5)... ',ratio);
        else
            options.Compression = false;
            print('skipped (compression ratio %.6g >= 0.5).\n',ratio);
        end
    else
        if options.Compression
            print('is mlsvd (requested by user)... ');
        else
            print('skipped (requested by user).\n');
        end
    end
elseif R == 1
    options.Compression = false;
    print('skipped (R = 1).\n');
else
    options.Compression = false;
    print('skipped (initialization U0 supplied).\n');
end
if options.Compression
    Tunc = T;
    [V,T] = lmlra_aca(T,size_core);
    if sum(size(T) > 1) >= 3
        output.Compression.ratio = ratio;
        output.Compression.relerr = frob(lmlrares(Tunc,V,T))/frob(Tunc);
        output.Compression.Name = 'Adaptive cross approximation';
        print('relative error = %.6g.\n',output.Compression.relerr);
    else
        options.Compression = false;
        T = Tunc; clear V;
        print('failed.\n');
    end
else
    output.Compression.Name = 'none';
end

% STEP 2: initialize the factor matrices unless they were all provided by
% the user.
print('Step 2: Initialization ');
if ~exist('U','var'), U = cell(1,N); end
Uempty = cellfun(@isempty,U);
if any(Uempty)
    if ischar(options.Initialization) && ...
       strcmpi(options.Initialization,'auto')
        if isnumeric(T) && sum(R <= size_tens) >= 2 && N > 2
            options.Initialization = @cpd_gevd;
            print('is cpd_gevd (sum(R <= size(T)) >= 2)... ');
        else
            options.Initialization = @cpd_rnd;
            print('is cpd_rnd (default)... ');
        end
    elseif isfunc(options.Initialization)
        print('is %s (requested by user)... ', ...
            func2str(options.Initialization));
    end
    if ~xsfunc(options.Initialization)
        error('cpd:Initialization','Not a valid initialization.');
    end
else
    options.Initialization = false;
    print('is manual... ');
end
if isfunc(options.Initialization)
    % Set the order in which the factor matrices should be updated.
    if sum(Uempty) > 0 && sum(Uempty) < N
        perm = [find(Uempty == 1) find(Uempty == 0)];
        if ~isfield(options,'AlgorithmOptions')
            options.AlgorithmOptions = struct;
        end
        if ~isfield(options.AlgorithmOptions,'Order')
            options.AlgorithmOptions.Order = perm;
        end
    end
    % Generate initial factor matrices.
    if ~isfield(options,'InitializationOptions')
        options.InitializationOptions = struct;
    end
    [U0,output.Initialization] = options.Initialization(T,R,...
                                 options.InitializationOptions);
    % Fill empty factor matrices.
    for n = 1:sum(Uempty), U{n} = U0{n}; end
    output.Initialization.Name = func2str(options.Initialization);
else
    output.Initialization.Name = 'manual';
end
output.Initialization.relerr = frob(cpdres(T,U))/frob(T);
print('relative error = %.6g.\n',output.Initialization.relerr);

% STEP 3: run the selected CPD algorithm.
if xsfunc(options.Algorithm)
    print('Step 3: Algorithm is %s... ',func2str(options.Algorithm));
    if ~isfield(options,'AlgorithmOptions')
        options.AlgorithmOptions = struct;
    end
    [U,output.Algorithm] = options.Algorithm(T,U,...
                           options.AlgorithmOptions);
    output.Algorithm.Name = func2str(options.Algorithm);
else
    error('cpd:Algorithm','Not a valid algorithm.');
end
if isfield(output.Algorithm,'iterations')
    print('iterations = %i, ',output.Algorithm.iterations);
end
output.Algorithm.relerr = frob(cpdres(T,U))/frob(T);
print('relative error = %.6g.\n',output.Algorithm.relerr);

% STEP 3b: expand to the original space.
if exist('V','var') && ~isempty(V)
    U = cellfun(@(u,v)v*u,U,V,'UniformOutput',false);
    T = Tunc;
end

% STEP 4: iteratively refine the solution if the tensor was compressed.
print('Step 4: Refinement ');
if ~options.Compression
    output.Refinement.Name = 'none';
    print('skipped (no decompression).\n');
elseif islogical(options.Refinement) && ~options.Refinement
    output.Refinement.Name = 'none';
    print('skipped (requested by user).\n');
elseif xsfunc(options.Refinement)
    print('is %s... ',func2str(options.Refinement));
    % Copy parameters from options.AlgorithmOptions which are not defined
    % in options.RefinementOptions.
    if ~isfield(options,'RefinementOptions')
        options.RefinementOptions = struct;
    end
    fn = fieldnames(options.AlgorithmOptions);
    for f = 1:length(fn)
        if ~isfield(options.RefinementOptions,fn{f})
            options.RefinementOptions.(fn{f}) = ...
                options.AlgorithmOptions.(fn{f});
        end
    end
    [U,output.Refinement] = options.Refinement(Tunc,U, ...
                            options.RefinementOptions);
    output.Refinement.Name = func2str(options.Refinement);
    if isfield(output.Refinement,'iterations')
        print('iterations = %i, ',output.Refinement.iterations);
    end
    output.Refinement.relerr = frob(cpdres(T,U))/frob(T);
    print('relative error = %.6g.\n',output.Refinement.relerr);
else
    error('cpd:Refinement','Not a valid refinement algorithm.');
end

% Format output.
fn = fieldnames(output);
for f = 1:length(fn)
    output.(fn{f}) = orderfields(output.(fn{f}));
end
