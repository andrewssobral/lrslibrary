function [U S V hist] = scaled_grass_mc(K, A_K, r, varargin)
% function [U S V hist] = scaled_grass_mc(K, A_K, k, varargin)
%
% Adapted from Boumal and Absil, 2011
% Adapted from C. T. Kelley, 1996-1997
%
% Matrix completion by geometric gradient descent with Armijo rule
% and polynomial linesearch 
% 
% Input: 
%   - K: *sparse* mask matrix of known entries (m x n)
%   - A_K: vector of known entries (nnz x 1)
%   - r: desired rank
%   - vargargin:
%     'maxit': max number of iteration (default 100)
%     'tol': tolerance (default 1.e-4)
%     'use_diag_s': diagonalized S during the iterations or not (default 0)
%     'beta_type': conjugate direction or not (default 'P-R')
%       - 'P-R': Polak Ribierre
%       - 'F-R': Fletcher Reeves
%       - 'Steep': Steepest descent
%     'grad_type': method to compute the gradients (default 'scaled')
%       - 'canon': canonical gradient
%       - 'scaled': scaled gradient
%       - 'sqrtscaled': scaled gradient with scaling = sqrt(SS')
%     'sigma_type': how to update S (default 'approx')
%       - 'approx': compute S approximately
%       - 'best': compute optimal S
%     'regval': regularize F (default 0.0)
%     'verbose': print out intermediate results? (default 0)
%     'rho': regularization term to keep incoherence (default 0.0)
%            (for most experiment just keep rho=0.0)
%     'U','V','S': Initializations of U, V, S
%            U,V,S will be initialized by svds(A) if these are not provided
% Output:
%   - U S V: U*S*V' is the approximation of A (U: mxr, S:rxr, V:nxr)
%   - hist: history during the iterations
%
% Thanh Ngo, 2012

params     = get_varargin(varargin);
maxit      = get_params(params, 'maxit', 100);
tol        = get_params(params, 'tol', 1.e-4);
use_diag_s = get_params(params, 'use_diag_s', 0);
beta_type  = get_params(params, 'beta_type', 'P-R');     % 'P-R', 'F-R' or 'Steep'
grad_type  = get_params(params, 'grad_type', 'scaled');  % 'canon', 'scaled', 'sqrtscaled'
sigma_type = get_params(params, 'sigma_type', 'approx'); % 'best', 'approx'
verbose    = get_params(params, 'verbose', 0);
%orth_value = get_params(params, 'orth_value', 0.1);
rho        = get_params(params, 'rho', 0.0);
regval     = get_params(params, 'regval', 0.0);

% make A_K a row vector
if size(A_K,1)>1
  A_K = A_K';
end

num_known = length(A_K);
[m n] = size(K);
[I J] = find(K);
MASK = sparse(I,J,ones(num_known,1),m,n);
I = int32(I);
J = int32(J);

% initialize U and V
if ~isfield(params,'U')
  % initialize with svds(A_K,0)
  [x.U x.S x.V] = svds(sparse(double(I),double(J),A_K,m,n),r);
else
  % use given initials U0 and V0
  x.U = params.U;
  x.V = params.V;
end

% compute S
if isfield(params,'S')
  x.S = params.S;
else
  x.S = getoptS_mex(MASK,x.U,x.V,A_K,I,J,regval);
end

% diagonalize S
if use_diag_s
  [R1 x.S R2] = svd(x.S);
  x.U = U*R1;
  x.V = V*R2;
end

start_time = tic;

x.US = x.U*x.S;
x.X = maskmult(x.US',x.V',I,J);
x.E = x.X-A_K;
x.Emat = sparse(double(I),double(J),x.E,m,n);

itc = 1;
xc = x; 
fc = F(xc,regval);
gc = grad(xc,rho,grad_type,0);
ip_gc = ip(xc,gc,gc,grad_type);
dir = scaleTxM(gc,-1); % first search-dir is steepest gradient

ithist= zeros(maxit,2);
step_sizes = zeros(maxit,1);
beta = 0;
total_time = 0;
for itc=1:maxit

  if itc>1
    start_time = tic;
  end

  % exact line search and back-tracking Armijo search
  step_sizes(itc) = exact_search_onR1(I,J,xc,dir);
  [xc_new,fc,succ,numf_a1,iarm1,step_sizes(itc)] = armijo_search(A_K,MASK,I,J,xc,fc,dir,step_sizes(itc),grad_type,sigma_type,regval);
  %[xc_new fc step_sizes(itc)] = backtrack(A_K,MASK,I,J,xc,fc,dir,sigma_type,regval);

  % gradients(x_new)
  gc_new = grad(xc_new,rho,grad_type,itc);
  ip_gc_new = ip(xc_new,gc_new,gc_new,grad_type);
  
  % test for convergence
  if sqrt(2*fc/num_known) < tol
    ithist(itc,1) = sqrt(2*fc/num_known); 
    total_time = total_time + toc(start_time);
    ithist(itc,2) = total_time;
    if verbose
      if isfield(params,'test_data')
        mae = mean(abs(maskmult(xc_new.S'*xc_new.U', xc_new.V', params.test_data.I, params.test_data.J)-params.test_data.T'));
        fprintf('Iter %d - Train RMSE %f - Test MAE: %f\n', itc, ithist(itc,1), mae);
      else
        fprintf('Iter %d - RMSE: %e\n', itc, ithist(itc,1));
      end
    end
    break;
  end

  % transport old gradient and search direction
  gc_old = transpVect(A_K,MASK,I,J,xc,gc,xc_new,1,grad_type);
  dir = transpVect(A_K,MASK,I,J,xc,dir,xc_new,1,grad_type);
  
  % we test how orthogonal the previous gradient is with
  % the current gradient, and possibly reset the to the gradient
  %orth_grads = ip(xc_new,gc_old,gc_new,grad_type)/ip_gc_new;
  
  % compute beta
  if strcmp(beta_type, 'steep') %| orth_grads >= orth_value
    dir = plusTxM(gc_new, dir, -1, 0);
  elseif strcmp(beta_type, 'F-R')  % Fletcher-Reeves
    beta = ip_gc_new/ip_gc
    dir = plusTxM(gc_new, dir, -1, beta);
  elseif strcmp(beta_type, 'P-R')  % Polak-Ribiere
    diff = plusTxM(gc_new, gc_old, 1, -1); % grad(new) - transported(grad(current))
    ip_diff = ip(xc_new, gc_new, diff, grad_type);
    beta = ip_diff / ip_gc;
    beta = max(0,beta);
    dir = plusTxM(gc_new, dir, -1, beta);
  end

  % update new to current
  gc = gc_new;
  ip_gc = ip_gc_new;
  xc = xc_new;

  ithist(itc,1) = sqrt(2*fc/num_known); 
  total_time = total_time + toc(start_time);
  ithist(itc,2) = total_time;
  if verbose
    if isfield(params,'test_data')
      mae = mean(abs(maskmult(xc_new.S'*xc_new.U', xc_new.V', params.test_data.I, params.test_data.J)-params.test_data.T'));
      fprintf('Iter %d - Train RMSE %f - Test MAE: %f\n', itc, ithist(itc,1), mae);
    else
      fprintf('Iter %d - RMSE: %e\n', itc, ithist(itc,1));
    end
  end

  if itc>5 & isfield(params,'tol_reschg')
    if abs(1-sqrt(ithist(itc,1)/ithist(itc-1,1)))<params.tol_reschg
      break;
    end
  end

end

U = xc_new.U;
S = xc_new.S;
V = xc_new.V;
hist.rmse = ithist(1:itc,1);
hist.time = ithist(1:itc,2);

%%---------- compute function F
function out = F(x, regval)
out = 0.5*(norm(x.E)^2);
if regval>0.0
  out = out + regval*(norm(x.S,'fro')^2);
end

%%---------- compute gradients
function g = grad(x, rho, grad_type, itc)
switch grad_type
case 'canon'
  EV = x.Emat*x.V*x.S';
  EtU = (x.US'*x.Emat)';
  g.U = EV - x.U*(x.U'*EV);
  g.V = EtU - x.V*(x.V'*EtU);
case 'scaled'
  if 0 & itc>10
    [U S V] = svd(x.S);
    s = diag(S); s = s.*(s>=1e-3) + 1e-3*(s<1e-3);
    IS = V*diag(1./s)*U';
  else
    IS = pinv(x.S);
  end
  EV = x.Emat*x.V*IS;
  EtU = (IS*x.U'*x.Emat)';
  g.U = EV - x.U*(x.U'*EV);
  g.V = EtU - x.V*(x.V'*EtU);
case 'sqrtscaled'
  EV = x.Emat*x.V;
  EtU = (x.U'*x.Emat)';
  g.U = EV - x.U*(x.U'*EV);
  g.V = EtU - x.V*(x.V'*EtU);
case 'scaled1'
  [U S V] = svd(x.S);
  temp = 1./diag(S); IS = diag(temp); IS2 = diag(temp.^2);
  D = V*IS*U'; D1 = U*IS2*U'; D2 = V*IS2*V';
  EV = x.Emat*x.V*D;
  EtU = (D*x.U'*x.Emat)';
  g.U = EV - x.U*(D1*(x.U'*EV));
  g.V = EtU - x.V*(D2*(x.V'*EtU));
case 'sqrtscaled1'
  [U S V] = svd(x.S);
  temp = 1./diag(S); IS = diag(temp);
  D1 = U*IS*U'; D2 = V*IS*V';
  EV = x.Emat*x.V;
  EtU = (x.U'*x.Emat)';
  g.U = EV - x.U*(D1*(x.U'*EV));
  g.V = EtU - x.V*(D2*(x.V'*EtU));
end
% regularization to keep incoherence
if rho>0.0
  fprintf('regularization');
  m0 = 10000;
  g.U = g.U + rho*Gp(x.U,m0,size(x.U,2));
  g.V = g.V + rho*Gp(x.V,m0,size(x.V,2));
end

%%---------- compute inner product
function i = ip(x,h1,h2,grad_type)
switch grad_type
case {'canon','scaled','sqrtscaled'}
  i = trace(h1.U'*h2.U) + trace(h1.V'*h2.V);
case 'scaled1'
  i = trace(x.S*x.S'*(h1.U'*h2.U)) + trace(x.S'*x.S*(h1.V'*h2.V));
case 'sqrtscaled1'
  [U S V] = svd(x.S);
  i = trace(U*S*U'*(h1.U'*h2.U)) + trace(V*S*V'*(h1.V'*h2.V));
end

%%---------- compute inner product
function h = scaleTxM(h,a)
h.U = a*h.U;
h.V = a*h.V;

function [xt ft t] = backtrack(A_K,K,I,J,xc,fc,gc,sigma_type,regval)
normg2 = norm(gc.U,'fro')^2 + norm(gc.V,'fro')^2;
normS2 = norm(xc.S,'fro')^2;
t = 1.e+2;
while t>1.e-8
  xt.U = xc.U + t*gc.U;
  xt.V = xc.V + t*gc.V;
  ft = 0.5*(norm(maskmult(xc.S'*xt.U',xt.V',I,J)-A_K)^2) + regval*normS2;
  if ft-fc < -.5*t*normg2
    break;
  end
  t = t/2;
end
xt = moveEIG(A_K,K,I,J,xc,gc,t,sigma_type,0,regval);
ft = F(xt,regval);


function tmin = exact_search_onR1(I,J,x,dir)
% Exact line search in the direction of dir on the tangent space of x
% !! so NOT the retracted curve, only use as guess !!
% x, current point
% dir is search direction
% returns: tmin

f0_omega = x.E;
dirUxSt = (dir.U*x.S)';
f1_omega = maskmult(dirUxSt,x.V',I,J)+maskmult(x.US',dir.V',I,J);
f2_omega = maskmult(dirUxSt,dir.V',I,J);

%a0 = norm(f0_omega,'fro')^2;
a1 = 2*sum(f0_omega.*f1_omega);
a2 = norm(f1_omega)^2 + 2*sum(f2_omega.*f0_omega);
a3 = 2*sum(f2_omega.*f1_omega);
a4 = norm(f2_omega)^2;

ts = roots([4*a4 3*a3 2*a2 a1]);
ind_real = find(ts==real(ts));
ts = ts(ind_real);
ind_real = find(ts>0);
if length(ind_real)>1
  tmin = max(ts(ind_real));
elseif length(ind_real)==0
  disp('no positive root');
  tmin = 1.0; % TODO
else
  tmin = ts(ind_real);
end

function [xt,ft,succ,numf,iarm,lambda] = armijo_search(A_K,K,I,J,xc,fc,gc,lambda,grad_type,sigma_type,regval)
% Armijo line search with polynomial interpolation
% Adapted from Boumal and Absil, 2011
% Adapted from steep by C. T. Kelley, Dec 20, 1996
% xc, current point
% fc = F(sys,xc)
% gc is search direction
%
% returns: xt,ft: Armijo point and F(sys,xt)
%          succ = true if success
%          numf: nb of F evals
%          iarm: nb of Armijo backtrackings

alp = 1e-4;
bhigh=.5; blow=.1;
MAX_IARM = 5;
numf = 0;

% trial evaluation at full step
xt = moveEIG(A_K,K,I,J,xc,gc,lambda,sigma_type,0,regval);
ft = F(xt,regval);
numf = numf+1; iarm = 0;

fgoal = fc-alp*lambda*ip(xc,gc,gc,grad_type);
% polynomial line search
q0=fc; qp0=-ip(xc,gc,gc,grad_type); lamc=lambda; qc=ft;
while(ft > fgoal)    
  iarm=iarm+1;
  if iarm==1
    lambda=polymod(q0, qp0, lamc, qc, blow, bhigh);
  else
    lambda=polymod(q0, qp0, lamc, qc, blow, bhigh, lamm, qm);
  end
  qm=qc; lamm=lamc; lamc=lambda;
  xt = moveEIG(A_K,K,I,J,xc,gc,lambda,sigma_type,0,regval);
  ft = F(xt,regval);
  numf = numf+1; qc=ft;
  if(iarm > MAX_IARM)
    succ = false;
    return
  end
  fgoal = fc-alp*lambda*ip(xc,gc,gc,grad_type);
end
succ = true;


function z = moveEIG(A_K,K,I,J,x,h,t,sigma_type,is_dir,regval)
%MOVEEIG  Retract a point x+t*h (minimal distance)
%  z = moveEIG(sys,x,h,t) retracts x+t*h to the manifold where
%     sys is a system,
%     x is a point of sys,
%     h is a tangent vector of TxM,
%     t is a distance (scalar).
[m n] = size(K);
[z.U Ru] = qr(x.U+t*h.U,0);
[z.V Rv] = qr(x.V+t*h.V,0);
if is_dir
  return;
end
switch sigma_type
case 'best'
  z.S = getoptS_mex(K,z.U,z.V,A_K,I,J,regval);
case 'approx'
  z.S = (z.U'*x.U)*x.S*(x.V'*z.V) - t*(z.U'*x.Emat*z.V);
end
z.US = z.U*z.S;
z.X = maskmult(z.US',z.V',I,J);
z.E = z.X-A_K;
z.Emat = sparse(double(I),double(J),z.E,m,n);

function [lplus]=polymod(q0, qp0, lamc, qc, blow, bhigh, lamm, qm)
%
% C. T. Kelley, Dec 29, 1997
%
% This code comes with no guarantee or warranty of any kind.
%
% function [lambda]=polymod(q0, qp0, qc, blow, bhigh, qm)
%
% Cubic/quadratic polynomial linesearch
%
% Finds minimizer lambda of the cubic polynomial q on the interval
% [blow * lamc, bhigh * lamc] such that
%
% q(0) = q0, q'(0) = qp0, q(lamc) = qc, q(lamm) = qm
% 
% if data for a cubic is not available (first stepsize reduction) then
% q is the quadratic such that
% 
% q(0) = q0, q'(0) = qp0, q(lamc) = qc
%
lleft=lamc*blow; lright=lamc*bhigh; 
if nargin == 6
%
% quadratic model (temp hedge in case lamc is not 1)
%
    lplus = - qp0/(2 * lamc*(qc - q0 - qp0) );
    if lplus < lleft lplus = lleft; end
    if lplus > lright lplus = lright; end
else
%
% cubic model
%
    a=[lamc^2, lamc^3; lamm^2, lamm^3];
    b=[qc; qm]-[q0 + qp0*lamc; q0 + qp0*lamm];
    c=a\b;
    lplus=(-c(1)+sqrt(c(1)*c(1) - 3 *c(2) *qp0))/(3*c(2));
    if lplus < lleft lplus = lleft; end
    if lplus > lright lplus = lright; end
end

function ht = transpVect(A_K,K,I,J,x,h,d,type,grad_type)
%TRANSPVECT  transport the vector (x,h) along the direction of d
% type == 0: transpVect(x,h,d)
%    we transport by projecting h orth. on R_x(d)
% type == 1: transpVect(x,h,z)
%    we transport by projecting h orth. on z
% returns ht, the transported h

%% z is the foot of ht
if nargin == 4
  type = 0;
end

if type==0    
  z = moveEIG(A_K,K,I,J,x,d,1,1,0.0);
else % type =1
  z = d;
end

%grad_type = 'canon';

%% the given vector is (x,h) and is transported onto d
ip_old = ip(x,h,h,grad_type);

switch grad_type
case {'canon','scaled','sqrtscaled'}
  ht.U = h.U - z.U*(z.U'*h.U);
  ht.V = h.V - z.V*(z.V'*h.V);
case 'scaled1'
  [U S V] = svd(x.S); temp = 1./diag(S); IS2 = diag(temp.^2);
  D1 = U*IS2*U';
  D2 = V*IS2*V';
  ht.U = h.U - z.U*(D1*(z.U'*h.U));
  ht.V = h.V - z.V*(D2*(z.V'*h.V));
case 'sqrtscaled1'
  [U S V] = svd(x.S); temp = 1./diag(S); IS = diag(temp);
  D1 = U*IS*U';
  D2 = V*IS*V';
  ht.U = h.U - z.U*(D1*(z.U'*h.U));
  ht.V = h.V - z.V*(D2*(z.V'*h.V));
end

ip_new = ip(z,ht,ht,grad_type);
ht = scaleTxM(ht, ip_old/ip_new);

function h = plusTxM(h1,h2,a1,a2)
% PLUSTXM  add 2 tangent vectors in the same tangent space
%  h = a1*h1 + a2*h2  a1,a2 reals
h.U = a1*h1.U + a2*h2.U;
h.V = a1*h1.V + a2*h2.V;

function out = Gp(X,m0,r)
z = sum(X.^2,2) /(2*m0*r) ;
z = 2*exp( (z-1).^2 ).*(z-1) ;
z(find(z<0)) = 0;
out = X.*repmat(z,1,r) / (m0*r) ;

