function [U,output] = cpd_als(T,U0,options)
%CPD_ALS CPD by alternating least squares.
%   [U,output] = cpd_als(T,U0) computes the factor matrices U{1}, ..., U{N}
%   belonging to a canonical polyadic decomposition of the N-th order
%   tensor T. The algorithm is initialized with the factor matrices U0{n}.
%   The structure output returns additional information:
%
%      output.alpha      - The value of the line or plane search step
%                          length(s) in every iteration.
%      output.fval       - The value of the objective function
%                          0.5*frob(T-cpdgen(U))^2 in every iteration.
%      output.info       - The circumstances under which the procedure
%                          terminated:
%                             1: Objective function tolerance reached.
%                             2: Step size tolerance reached.
%                             3: Maximum number of iterations reached.
%      output.iterations - The number of iterations.
%      output.relfval    - The difference in objective function value
%                          between every two successive iterates, relative
%                          to its initial value.
%      output.relstep    - The step size relative to the norm of the 
%                          current iterate in every iteration.
%
%   cpd_als(T,U0,options) may be used to set the following options:
%
%      options.Display = 0        - Displays output information each
%                                   options.Display iterations.
%      options.LineSearch =       - A function handle to the desired line
%      [{false}|@cpd_aels|...       search algorithm.
%       @cpd_els|@cpd_lsb]
%      options.LineSearchOptions  - An options structure passed to the
%                                   selected line search algorithm.
%      options.PlaneSearch =      - A function handle to the desired plane
%      [{false}|@cpd_eps]           search algorithm. Searches in the plane
%                                   spanned by the last two ALS updates.
%      options.PlaneSearchOptions - An options structure passed to the
%                                   selected plane search algorithm.
%      options.MaxIter = 500      - The maximum number of iterations.
%      options.Order = 1:N        - The order in which to update the factor
%                                   matrices. Must be a permutation of 1:N.
%      options.TolFun = 1e-8      - The tolerance for output.relfval. Note
%                                   that because the objective function is
%                                   a squared norm, TolFun can be as small
%                                   as eps^2.
%      options.TolX = 1e-6        - The tolerance for output.relstep.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Format the tensor T.
T = fmt(T,true);

% Check the initial factor matrices U0.
U = U0(:).'; U1 = {};
R = size(U{1},2);
N = length(U);
if any(cellfun('size',U,2) ~= R)
    error('cpd_als:U0','size(U0{n},2) should be the same for all n.');
end

% Check the options structure.
isfunc = @(f)isa(f,'function_handle');
xsfunc = @(f)isfunc(f)&&exist(func2str(f),'file');
if nargin < 3, options = struct; end
if ~isfield(options,'Display'), options.Display = 0; end
if ~isfield(options,'FastUpdate'), options.FastUpdate = true; end
if ~isfield(options,'LineSearch'), options.LineSearch = false; end
if ~xsfunc(options.LineSearch) ...
   && (~isa(options.LineSearch,'function_handle') && options.LineSearch)
    error('cpd_als:LineSearch','Not a valid line search algorithm.');
end
if ~isfield(options,'LineSearchOptions')
    options.LineSearchOptions = struct;
end
if ~isfield(options,'PlaneSearch'), options.PlaneSearch = false; end
if ~xsfunc(options.PlaneSearch) ...
   && (~isa(options.PlaneSearch,'function_handle') && options.PlaneSearch)
    error('cpd_als:PlaneSearch','Not a valid line search algorithm.');
end
if ~isfield(options,'PlaneSearchOptions')
    options.PlaneSearchOptions = struct;
end
if ~isfield(options,'MaxIter'), options.MaxIter = 500; end
if ~isfield(options,'Order'), options.Order = 1:N; end
if ~isfield(options,'TolFun'), options.TolFun = 1e-8; end
if ~isfield(options,'TolX'), options.TolX = 1e-6; end

% Cache some intermediate variables.
T2 = frob(T,'squared');
K = cell(1,N);
UHU = zeros(R,R,N);
for n = 1:N, UHU(:,:,n) = U{n}'*U{n}; end

% Alternating least squares.
first = options.Order(1);
last = options.Order(end);
K{first} = mtkrprod(T,U,first);
D = cpdres(T,U);
output.alpha = [];
output.fval = 0.5*(D(:)'*D(:));
output.info = false;
output.iterations = 0;
output.relgain = [];
output.relfval = [];
output.relstep = [];
while ~output.info
    
    % Save current state.
    U2 = U1;
    U1 = U;
    
    % Update factor matrices.
    for n = options.Order
        W = prod(UHU(:,:,[1:n-1 n+1:N]),3);
        if n ~= first
            K{n} = mtkrprod(T,U,n);
        end
        U{n} = K{n}/conj(W);
        UHU(:,:,n) = U{n}'*U{n};
    end
    
    % Line/plane search.
    [U,alpha,outputsrch] = search();
	
    % Update the output structure.
    output.alpha(:,end+1) = alpha;
    if isfield(outputsrch,'relgain')
        output.relgain(end+1) = outputsrch.relgain;
    end
    output.fval(end+1) = outputsrch.fval;
    output.iterations = output.iterations+1;
    output.relfval(end+1) = ...
        abs(diff(output.fval(end:-1:end-1)))/abs(output.fval(1));
    output.relstep(end+1) = sqrt(sum(cellfun(@(u,v)(u(:)-v(:))'* ...
        (u(:)-v(:)),U,U1)))/sqrt(sum(cellfun(@(u)u(:)'*u(:),U)));
    if output.relfval(end) <= options.TolFun, output.info = 1; end
    if output.relstep(end) <= options.TolX, output.info = 2; end
    if output.iterations >= options.MaxIter, output.info = 3; end
    
    % Display progress.
    if options.Display > 0 && (output.iterations == 1 || output.info || ...
       mod(output.iterations,options.Display) == 0)
        if output.iterations == 1
            bold = '%s';
            [~,~,~,~,v] = regexp(version('-release'),'([0-9]+)([ab])');
            if usejava('Desktop') && str2double(v{1}{1}) > 2011 || ...
               (str2double(v{1}{1}) == 2011 && strcmpi(v{1}{2},'b'))
                bold = '<strong>%s</strong>';
            end
        end
        if output.iterations == 1 || ...
           mod(output.iterations,15*options.Display) == 0
            fprintf('\n%7s%s','',sprintf(bold,'fval'));
            fprintf('%13s%s','',sprintf(bold,'relfval'));
            fprintf('%10s%s','',sprintf(bold,'relstep'));
            fprintf('%10s%s','',sprintf(bold,'alpha'));
            fprintf('\n%30s = %4.e %6s = %4.e\n\n', ...
                    'TolFun',options.TolFun,'TolX',options.TolX);
        end
        if output.iterations == 1
            fprintf('%4i: % 14.8e |\n',0,output.fval(1));
        end
        stralpha = [repmat('%10.4e ',1,size(output.alpha,1)) '\n'];
        fprintf(['%4i: % 14.8e | %14.8e | %14.8e | ' stralpha], ...
                output.iterations,output.fval(end), ...
                output.relfval(end),output.relstep(end), ...
                abs(output.alpha(:,end)));
    end

end

% Display termination message.
if options.Display > 0
    ahref = '\n%s\n\n';
    x = round(linspace(0,output.iterations,min(500,output.iterations)));
    if length(bold) > 2
        ahref = sprintf(['\n<a href="matlab:semilogy(%s,%s);' ...
            'xlabel(''iteration'');legend(''fval'',' ...
            '''relfval'',''relstep'')">%%s</a>\n\n'],mat2str(x'), ...
            mat2str([output.fval(x+1)' [nan output.relfval(x(2:end))]' ...
                    [nan output.relstep(x(2:end))]'],3));
    end
    switch output.info
        case 1, fprintf(ahref,'Objective function tolerance reached.');
        case 2, fprintf(ahref,'Step size tolerance reached.');
        case 3, fprintf(ahref,'Maximum number of iterations reached.');
    end
end

function [Ua,alpha,outputsrch,dU,dU1,dU2,D,n] = search()
% Line/plane search wrapper that guarantees the return values alpha and
% outputsrch.fval, and prepares the next ALS iteration.
    
    % Attempt line/plane search, and fill in missing alpha.
    outputsrch = struct;
    state = output; state.T2 = T2; state.UHU = UHU;
    if options.FastUpdate
        state.K1 = first; state.KN = last; state.K = K;
    end
    if isfunc(options.LineSearch)
        dU = cellfun(@(u,v)u-v,U,U1,'UniformOutput',false);
        [alpha,outputsrch] = options.LineSearch(T,U1,dU, ...
            state,options.LineSearchOptions);
        if any(isnan(alpha)) || isempty(alpha), alpha = [1 1];
        elseif length(alpha) == 1, alpha(2) = 1; end
    elseif isfunc(options.PlaneSearch)
        if output.iterations >= 1
            dU1 = cellfun(@(u,v)u-v,U,U1,'UniformOutput',false);
            dU2 = cellfun(@(u,v)u-v,U1,U2,'UniformOutput',false);
            [alpha,outputsrch] = options.PlaneSearch(T,U1,dU1,dU2, ...
                state,options.PlaneSearchOptions);
        else
            alpha = nan;
        end
        if any(isnan(alpha)) || isempty(alpha), alpha = [1 0 1];
        elseif length(alpha) == 2, alpha(3) = 1; end
    else
        alpha = 1;
    end
    
    % Fill in missing fval.
    if ~isfield(outputsrch,'fval')
        if isfunc(options.LineSearch)
            Ua = cellfun(@(u,v)alpha(2)*(u+alpha(1)*v),U1,dU, ...
                'UniformOutput',false);
        elseif isfunc(options.PlaneSearch) && output.iterations >= 1
            Ua = cellfun(@(u,v,w)alpha(3)*(u+alpha(1)*v+alpha(2)*w), ...
                U1,dU1,dU2,'UniformOutput',false);
        else
            Ua = U;
        end
        if options.FastUpdate
            K1 = mtkrprod(T,Ua,first);
            UHU1 = zeros(size(UHU));
            for n = 1:N, UHU1(:,:,n) = Ua{n}'*Ua{n}; end
        end
        if options.FastUpdate && log10(output.fval(end)) > log10(T2)-16+2.5
            outputsrch.fval = abs(.5*(T2+sum(sum(real(prod(UHU1,3)))))- ...
                            real(sum(dot(K1,Ua{first}))));
        else
            D = cpdres(T,Ua);
            outputsrch.fval = 0.5*(D(:)'*D(:));
        end
        if outputsrch.fval > output.fval(end), K1 = []; end
    end
    
    % If line search did not improve the objective function, don't take it.
    if outputsrch.fval > output.fval(end)
        if isfunc(options.LineSearch)
            alpha = [1 1];
        elseif isfunc(options.PlaneSearch) && output.iterations >= 1
            alpha = [1 0 1];
        else
            alpha = 1;
        end
        if options.FastUpdate && log10(output.fval(end)) > log10(T2)-16+2.5
            outputsrch.fval = abs(0.5*(T2+sum(sum(real( ...
                W.*UHU(:,:,last)))))-real(sum(dot(K{last},U{last}))));
        else
            D = cpdres(T,U);
            outputsrch.fval = 0.5*(D(:)'*D(:));
        end
    end
    
    % Prepare next ALS iteration, if necessary.
    if options.FastUpdate
        if exist('K1','var') && ~isempty(K1)
            K{first} = K1;
            UHU = UHU1;
        else
            if isfunc(options.LineSearch)
                Ua = cellfun(@(u,v)alpha(2)*(u+alpha(1)*v),U1,dU, ...
                    'UniformOutput',false);
            elseif isfunc(options.PlaneSearch) && output.iterations >= 1
                Ua = cellfun(@(u,v,w)alpha(3)*(u+alpha(1)*v+alpha(2)*w),...
                    U1,dU1,dU2,'UniformOutput',false);
            else
                Ua = U;
            end
            K{first} = mtkrprod(T,Ua,first);
            for n = 1:N, UHU(:,:,n) = Ua{n}'*Ua{n}; end
        end
    end
    
end

end
