function [A,Out] = ncp(M,r,opts)
% ncp: nonnegative tensor decomposition (CANDECOMP/PARAFAC) by block-coordinate update
%  min 0.5*||M - A_1\circ...\circ A_N||_F^2  
%  subject to A_1>=0, ..., A_N>=0
%
% input: 
%       M: input nonnegative tensor
%       r: estimated rank (each A_i has r columns); require exact or moderate overestimates
%       opts.
%           tol: tolerance for relative change of function value, default: 1e-4
%           maxit: max number of iterations, default: 500
%           maxT: max running time, default: 1e3
%           rw: control the extrapolation weight, default: 1
%           A0: initial point in cell struct, default: Gaussian random
%           matrices
% output:
%       A: nonnegative ktensor
%       Out.
%           iter: number of iterations
%           hist_obj: history of objective values
%           hist_rel: history of relative objective changes (row 1) and relative residuals (row 2)
%
% require MATLAB Tensor Toolbox from
% http://www.sandia.gov/~tgkolda/TensorToolbox/
%
% More information can be found at:
% http://www.caam.rice.edu/~optimization/bcu/

%% Parameters and defaults
if isfield(opts,'tol');   tol = opts.tol;     else tol = 1e-4;   end % stopping tolerance
if isfield(opts,'maxit'); maxit = opts.maxit; else maxit = 500;  end % max # of iterations
if isfield(opts,'maxT');  maxT = opts.maxT;   else maxT = 1e3;   end % max time in seconds 
if isfield(opts,'rw');    rw = opts.rw;       else rw = 1;       end % initial extrapolation weight

%% Data preprocessing and initialization
N = ndims(M); % M is an N-way tensor 
nWay = M.size; % dimensions of M
Mnrm = norm(M); % norm of M
obj0 = 0.5*Mnrm^2; % initial objective value

% initial tensor factors
if isfield(opts,'A0') 
    A0 = opts.A0; 
else
    A0 = cell(1,N);
    for n = 1:N
        A0{n} = max(0,randn(nWay(n),r)); % randomly generate each factor
    end
end

% normalize A0 and cache its square
Asq = cell(1,N);
for n = 1:N
    A0{n} = A0{n}/norm(A0{n},'fro')*Mnrm^(1/N);
    Asq{n} = A0{n}'*A0{n};
end
Am = A0; A = A0; 

nstall = 0; % # of stalled iterations
t0 = 1; % used for extrapolation weight update
wA = ones(N,1); % extrapolation weight array
L0 = ones(N,1); L = ones(N,1); % Lipschitz constant array

%% Store data in memory to save time if the data is not too large
if N*prod(nWay)<4000^2
    storedata = true;
    sizeN = zeros(1,N);
    pM = cell(1,N); kroA = cell(1,N);
    % store every mode-matricization of M
    for n = 1:N
        pM{n} = permute(M,[n 1:n-1,n+1:N]);
        sizeN(n) = prod(nWay)/size(M,n);
        pM{n} = reshape(pM{n}.data, size(M,n), sizeN(n));
        kroA{n} = zeros(sizeN(n),r);
    end
else
    storedata = false;
end

%% Iterations of block-coordinate update 
%
%  iteratively updated variables:
%       Gn: gradients with respect to A{n}
%       A: new updates
%       A0: old updates
%       Am: extrapolations of A
%       L, L0: current and previous Lipschitz bounds
%       obj, obj0: current and previous objective values

start_time = tic;
fprintf('Iteration:     ');

for k = 1:maxit
    fprintf('\b\b\b\b\b%5i',k);    
    
    %--- update each factor matrix ----%
    % For derivation of the update, see Section 3.2 of this report:
    % Yangyang Xu and Wotao Yin. A block coordinate descent method for 
    % regularized multi-convex optimization with application to nonnegative
    % tensor factorization and completion. 
    % Rice University CAAM Technical Report TR12-15
    for n = 1:N
        % Bsq = Asq{1}.*Asq{2}...Asq{n-1}.*Asq{n+1}....*Asq{N}
        Bsq = ones(r);
        for i = 1:N
            if i ~= n
                Bsq = Bsq.*Asq{i};
            end
        end
        L0(n) = L(n);
        L(n) = norm(Bsq); % gradient Lipschitz constant
        if storedata
            % do Khatri-Rao-product
            % kroA{n} = A{N}\odot...\odot A{n+1}\odot A{n-1}...\odot A{1}
            matorder = [N:-1:n+1,n-1:-1:1];
            for j = 1:r
                ab = A{matorder(1)}(:,j);
                for i = matorder(2:N-1)
                    ab = A{i}(:,j) * ab(:).';
                end
                kroA{n}(:,j) = ab(:);
            end
            % do matricized-tensor-times-Khatri-Rao-product
            % which equals n-th mode matricization of M times kroA
            MB = pM{n}*kroA{n};
        else
            % do matricized-tensor-times-Khatri-Rao-product
            MB = mttkrp(M,A,n);
        end
        %compute the gradient
        Gn = Am{n}*Bsq-MB;
        A{n} = max(0,Am{n}-Gn/L(n));
        Asq{n} = A{n}'*A{n};
    end
    
    % --- diagnostics, reporting, stopping checks ---
    obj = 0.5*(sum(sum(Asq{N}.*Bsq))-2*sum(sum(A{N}.*MB))+Mnrm^2); % current objective value
    relerr1 = abs(obj-obj0)/(obj0+1); % relative objective change
    relerr2 = (2*obj)^.5/Mnrm; % relative residual
    
    % reporting
    Out.hist_obj(k) = obj;
    Out.hist_rel(1,k) = relerr1;      
    Out.hist_rel(2,k) = relerr2;
    
    % check stopping criterion
    crit = relerr1<tol; 
    if crit; nstall = nstall+1; else nstall = 0; end
    if nstall>=3 || relerr2<tol; break; end
    if toc(start_time)>maxT; break; end;
    
    % --- correction and extrapolation ---
    t = (1+sqrt(1+4*t0^2))/2;
    if obj>=obj0 
        % restore previous A to make the objective nonincreasing
        Am = A0;
    else
        % apply extrapolation
        w = (t0-1)/t; % extrapolation weight
        for n = 1:N
            wA(n) = min([w,rw*sqrt(L0(n)/L(n))]); % choose smaller weights for convergence
            Am{n} = A{n}+wA(n)*(A{n}-A0{n}); % extrapolation
        end
        A0 = A; t0 = t; obj0 = obj;
    end    
end
fprintf('\n');
A = ktensor(A);
Out.iter = k; % report # of iterations

%% <http://www.caam.rice.edu/~optimization/bcu/ncp/ncp.m Download this m-file>