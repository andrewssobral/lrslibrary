function [A,C,Out] = ntd(M,coreNway,opts)
% ntd: nonnegative Tucker decomposition by block-coordinate
% update
%  min 0.5*||M - C \times_1 A_1 ...\times_N A_N||_F^2  
%  subject to A_1>=0, ..., A_N>=0, C>=0
% input: 
%       M: input nonnegative tensor
%       coreNway: size of core tensor C
%       opts.
%           tol: tolerance for relative change of function value, default:
%           1e-4
%           maxit: max number of iterations, default: 500
%           maxT: max running time, default: 1e3
%           rw: control the extrapolation weight, default: 1
%           A0: initial point in cell struct, default: Gaussian random
%           matrices
%           C0: initial value of C, default: Gaussian random array
% output:
%       A: cell struct with each component being nonnegative matrix
%       C: nonnegative core tensor
%       Out.
%           iter: number of iterations
%           hist_obj: history of objective values at each iteration
%           hist_rel: history of relative changes at each iteration
%
% require the Toolbox of tensor
% downloaded from http://www.sandia.gov/~tgkolda/TensorToolbox/
%
% More information can be found at:
% http://www.caam.rice.edu/~optimization/bcu/

%% Parameters and defaults
if isfield(opts,'maxit')  maxit = opts.maxit; else maxit = 500;  end
if isfield(opts,'tol')    tol = opts.tol;     else tol = 1e-4;   end
if isfield(opts,'rw')     rw = opts.rw;       else rw = 1;       end
if isfield(opts,'maxT')   maxT = opts.maxT;   else maxT = 1e3;   end

%% Data preprocessing and initialization
N = ndims(M); % M is an N-way tensor
Nway = M.size; % dimension of M

if isfield(opts,'A0')
    A0 = opts.A0;
else
    A0 = cell(1,N);
    for n = 1:N
        % randomly generate each factor matrix
        A0{n} = max(0,randn(Nway(n),coreNway(n)));
    end
end

if isfield(opts,'C0')  
    C0 = opts.C0; 
else
    % randomly generate core tensor
    C0 = tensor(max(0,randn(coreNway)));
end

Mnrm = norm(M); obj0 = 0.5*Mnrm^2;

Asq = cell(1,N);
for n = 1:N
    A0{n} = A0{n}/norm(A0{n},'fro')*Mnrm^(1/(N+1));
    Asq{n} = A0{n}'*A0{n};
end
C0 = tensor(C0/norm(C0)*Mnrm^(1/(N+1)));
A = A0; Am = A0;
C = C0; Cm = C0;

nstall = 0; w = 0; t0 = 1; t = t0; wA = ones(N+1,1);
L0 = ones(N+1,1); L = ones(N+1,1);

%% Store data to save computing time if it is not too large
storedata = false;
if N*prod(Nway)<4000^2
    storedata = true;
    pM = cell(1,N);
    for n = 1:N
        pM{n} = permute(M,[n 1:n-1,n+1:N]);
        sizeN = prod(Nway)/size(M,n);
        pM{n} = reshape(pM{n}.data, size(M,n), sizeN);
    end
end

%% Iterations of block-coordinate update 
%
%  iteratively updated variables:
%       GradA: gradients with respect to each component matrix of A
%       GradC: gradient with respect to C
%       A,C: new updates
%       A0,C0: old updates
%       Am,Cm: extrapolations of A
%       L, L0: current and previous Lipschitz bounds
%       obj, obj0: current and previous objective values

start_time = tic;
fprintf('Iteration:     ');

for k = 1:maxit
    fprintf('\b\b\b\b\b%5i',k);
    
    % -- update the core tensor C --
    L0(N+1) = L(N+1);
    L(N+1) = 1;
    for i = 1:N L(N+1) = L(N+1)*norm(Asq{i}); end
    Bsq = Cm; GradC = M;
    for i = 1:N
        Bsq = ttm(Bsq,Asq{i},i);
        GradC = ttm(GradC,A{i}',i);
    end
    %compute the gradient
    GradC = Bsq-GradC;
    C = tensor(max(0,Cm.data-GradC.data/L(N+1)));
    
    % -- update factor matrices A --
    for n = 1:N
        if storedata
            B = C;
            for i = 1:N
                if i~=n
                    B = ttm(B,A{i},i);
                end
            end
            B = tenmat(B,n);
            Bsq = B.data*B.data'; 
            MB = pM{n}*B.data';
        else
            Bsq = C;
            for i = 1:N
                if i~=n
                    Bsq = ttm(Bsq,Asq{i},i);
                end
            end
            Bsq = tenmat(Bsq,n)*tenmat(C,n)';          
            
            MB = M;
            for i = 1:N
                if i~=n
                    MB = ttm(MB,A{i}',i);
                end
            end
            MB = tenmat(MB,n)*tenmat(C,n)';
            Bsq = Bsq.data; MB = MB.data;
        end 
        %compute the gradient
        GradA = Am{n}*Bsq-MB; 
        L0(n) = L(n);
        L(n) = norm(Bsq);
        A{n} = max(0,Am{n}-GradA/L(n));
        Asq{n} = A{n}'*A{n};
    end
    
    % --- diagnostics, reporting, stopping checks ---
    obj = 0.5*(sum(sum(Asq{N}.*Bsq))-2*sum(sum(A{N}.*MB))+Mnrm^2);
    relerr1 = abs(obj-obj0)/(obj0+1);    relerr2 = (2*obj)^.5/Mnrm;
    
    % reporting
    Out.hist_obj(k) = obj;
    Out.hist_rel(1,k) = relerr1;      
    Out.hist_rel(2,k) = relerr2;
    
    % check stopping criterion   
    crit = relerr1<tol;
    if crit; nstall = nstall+1; else nstall = 0; end
    if nstall>=3 || relerr2<tol break; end
    if toc(start_time)>maxT; break; end;
    
    % --- correction and extrapolation ---
    t = (1+sqrt(1+4*t0^2))/2;
    if obj>=obj0
        % restore A, C to make the objective nonincreasing
        Am = A0; Cm = C0;
    else
        %get new extrapolated points
        w = (t0-1)/t; % extrapolation weight
        for n = 1:N
            % choose smaller weight for convergence
            wA(n) = min([w,rw*sqrt(L0(n)/L(n))]);
            Am{n} = A{n}+wA(n)*(A{n}-A0{n});
        end
        % choose smaller weight for convergence
        wA(N+1) = min([w,rw*sqrt(L0(N+1)/L(N+1))]);
        Cm = tensor(C.data+wA(N+1)*(C.data-C0.data));
        A0 = A; C0 = C; t0 = t; obj0 = obj;
    end
end
fprintf('\n'); Out.iter = k;

%% <http://www.caam.rice.edu/~optimization/bcu/ntd/ntd.m Download this m-file>