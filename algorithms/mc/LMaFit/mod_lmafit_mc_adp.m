function [X,Y, infos] = mod_lmafit_mc_adp(m,n,Known,data_ts,opts)
%
% Solver for matrix completion:
%
%        U_ij ~= A_ij,  for i,j in Known
%
% where U =  X*Y and A_ij's are given in an known index set.
%
% Output:
%           X --- m x k matrix
%           Y --- k x n matrix
%         Out --- output information
% Input:
%        m, n --- matrix sizes
%           k --- rank estimate
%        data --- values of known elements in a 1D row vector
%       Known --- positions of known elements in a 1D row vector
%                 assuming matrices are arranged column-wise
%             or  Known is a structure with two fields [Ik, Jk]
%        opts --- option structure with fields: opt, maxit
%
% Copyright(c) 2009 Yin Zhang
%
%   modified by Zaiwen Wen, 12/17/2009
%


% set parameters
t_begin = tic();

tol = 1.25e-4;
maxit = 500;

if isfield(opts,'tol');         tol     = opts.tol;        end
if isfield(opts,'maxit');       maxit   = opts.maxit;      end

I = Known.Ik; J = Known.Jk; data = Known.entries;
datanrm = max(1,norm(data));
S = sparse(I, J, data, m, n);
[I, J, data] = find(S);
data = data';


L = length(data); omega = 1; omega_max = 50; delta = 1; gamma = 0.7;


X = opts.X;
Y = opts.Y;
Res = data - partXY(X',Y,I,J,L);
%  updateSval(S, Res, L);
updateSparse(S, Res);
resnrm = norm(Res);
res = data;

infos.costs = (resnrm^2) /L;
infos.iter_time = 0;
% If interested in computing recovery
if opts.compute_predictions,
    preds_test = partXY(X', Y, data_ts.rows,data_ts.cols, data_ts.nentries)';
    errors_test = preds_test - data_ts.entries;
    cost_test = (errors_test'*errors_test)/data_ts.nentries;
    infos.test_error = cost_test;
end
if opts.verbosity,
    fprintf('iter: %3i cost: %3e \n',0, (resnrm^2) / L);
end;
% main iteration

for iter = 1:maxit
    %% step 4
    
    % keep the old values
    Xo = X; Yo = Y; res0 = res; resnrm0 = resnrm; omega0 = omega;
    
    % (2.11a) using the notation of (2.10)
    SY = S*Y';
    X = SY + X*(Y*Y');
    [X, ~] = qr(X,0);
    
    % (2.11b) using the notation of (2.10)
    XS = X'*S;
    Y = XS + (X'*Xo)*Y;
    
    % (2.11c) & (2.11d) using the notation of (2.10)
    dd = partXY(X',Y,I,J, L);
    res = data - dd;
    
    % compute the residual ratio
    resnrm = norm(res);  relres = resnrm/datanrm ;  ratio = resnrm/resnrm0;
    
    % for 'stagnation stopping criterion'
    reschg = abs(1-resnrm/resnrm0);
    %% step 6
    
    if ratio >= 1
        delta = max(0.1*(omega-1), 0.1*delta);
        X = Xo; Y = Yo; res = res0; resnrm = resnrm0; relres = resnrm/datanrm ;
        omega = 1;
        
        %% step 8
        
    elseif ratio > gamma
        delta = max(delta,0.25*(omega-1));
        omega = min(omega + delta,omega_max);
    end
    
    %     updateSval(S,omega*res,L);
    updateSparse(S,omega*res);
    obj_val = (resnrm^2) /L;
    
    
    % printout
    if opts.verbosity,
        fprintf('iter: %3i cost: %7.3e relres: %3.1e ratio: %4.4f chg: %3.1e omega: %3.1e delta: %3.1e\n',...
            iter,obj_val, relres,ratio,reschg,omega0,delta);
    end
    
    infos.costs = [infos.costs; obj_val];
    infos.iter_time = [infos.iter_time; toc(t_begin)];
    % If interested in computing recovery
    if opts.compute_predictions,
        preds_test = partXY(X', Y, data_ts.rows,data_ts.cols, data_ts.nentries)';
        errors_test = preds_test - data_ts.entries;
        cost_test = (errors_test'*errors_test)/data_ts.nentries;
        infos.test_error = [infos.test_error; cost_test];
    end
    if obj_val < tol,
        break;
    end
    
    
end %iter

infos.time = sum(infos.iter_time);

end %main
