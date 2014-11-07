function [U,S,output] = lmlra3_rtr(T,U0,options)
%LMLRA3_RTR LMLRA by a Riemannian trust region method.
%   [U,S,output] = lmlra3_rtr(T,U0) computes the factor matrices U{1}, 
%   U{2} and U{3} and core tensor S belonging to a low multilinear rank
%   approximation of the third order tensor T. The algorithm is initialized
%   with the factor matrices U0{n}. The structure output returns additional
%   information:
%
%      output.accepted - Whether the proposed iterate was accepted or not.
%      output.cauchy   - Whether the Cauchy point was used in place of an
%                        update computed from a random initial tangent
%                        vector.
%      output.Delta    - The trust-region radius at the outer iteration.
%      output.dist     - The distance from the solution.
%      output.normS    - The Frobenius norm of the core tensor S. The
%                        objective function is to minimize -frob(S)^2.
%      output.k        - The outer iteration number.
%      output.ng       - The norm of the gradient.
%      output.numinner - The number of inner iterations used to compute the
%                        next iterate.
%      output.rho      - The performance ratio for the iterate.
%      output.sangle   - The difference in subspace angle of U{1} between
%                        every two successive iterates.
%      output.time     - The wallclock time for the outer iteration.
%
%   lmlra3_rtr(T,U0,options) may be used to set the following options:
%
%      options.Delta_bar = inf - Maximum trust-region radius.
%      options.Delta0          - Initial trust-region radius.
%      = sum(size(U0,2))
%      options.epsilon = 1e-6  - Outer convergence tolerance (absolute).
%      options.kappa = 0.1     - Inner kappa convergence tolerance.
%      options.max_inner       - Maximum number of inner iterations.
%      = dot(size(T)-size(U0,2),size(U0,2))
%      options.max_outer = 100 - Maximum number of outer iterations.
%      options.min_inner = 0   - Minimum number of inner iterations.
%      options.min_outer = 0   - Minimum number of outer iterations.
%      options.rho_prime = 0.1 - Accept/reject ratio.
%      options.theta = 1.0     - Inner theta convergence tolerance.
%      options.verbosity = 0   - Level of verbosity:
%                                   0: Silent.
%                                   1: One-line info.
%                                   2: Detailed info.

%   Authors: Mariya Ishteva (mariya.ishteva@cc.gatech.edu)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   Based on a modified version of GenRTR (Generic Riemannian Trust Region
%   Package), copyright P.-A. Absil, C.G. Baker, and K.A. Gallivan, Florida
%   State University, School of Computational Science.
%
%   References:
%   [1] M. Ishteva, P.-A. Absil, S. Van Huffel, L. De Lathauwer, "Best low
%       multilinear rank approximation of higher-order tensors, based on
%       the Riemannian trust-region scheme", SIAM J. Matrix Anal. Appl.,
%       Vol. 32, No. 1, 2011, pp. 115-135.

% Check the tensor T.
N = ndims(T);
if N ~= 3, error('lmlra3_rtr:T','ndims(T) should be 3.'); end
if ~isreal(T), error('lmlra3_rtr:T','isreal(T) should be true.'); end

% Check the initial factor matrices U0.
if any(cellfun('size',U0(:).',1) ~= size(T))
    error('lmlra3_rtr:U0','size(T,n) should equal size(U0{n},1).');
end

% Set up dimensions and initial iterate.
options.x0.U = U0{1};
options.x0.V = U0{2};
options.x0.W = U0{3};
[I1,I2,I3] = size(T);
R1 = size(U0{1},2);
R2 = size(U0{2},2);
R3 = size(U0{3},2);
d = (I1-R1)*R1+(I2-R2)*R2+(I3-R3)*R3;

% Check the options structure.
if ~isfield(options,'max_inner'), options.max_inner = d; end
if options.max_inner > d, options.max_inner = d; end
if ~isfield(options,'Delta_bar'), options.Delta_bar = inf; end
if ~isfield(options,'Delta0'), options.Delta0 = R1+R2+R3; end
if ~isfield(options,'verbosity'), options.verbosity = 0; end

% Set up function handles.
fns.R     = @(x,eta)R(T,x,eta); % Retraction
fns.g     = @g;                 % Riemannian metric
fns.proj  = @proj;              % Projection onto tangent plane from nearby
fns.randT = @randT;
fns.f     = @(x)f(T,x);         % Objective function
fns.fgrad = @(x)grad(x);        % Gradient of f
fns.fhess = @(x,eta)H(x,eta);   % Hessian of f
Tzero.U = zeros(I1,R1);
Tzero.V = zeros(I2,R2);
Tzero.W = zeros(I3,R3);
options.x0 = R(T,options.x0,Tzero);
clear Tzero;

% Call Riemannian trust region solver.
[x,stats] = rtr_Gr3(fns,options);
U = {x.U,x.V,x.W};
S = x.AUVW;
fn = fieldnames(stats);
for i = 1:length(fn)
    output.(fn{i}) = [];
    for n = 1:length(stats)
        output.(fn{i}) = [output.(fn{i}) stats(n).(fn{i})];
    end
end

end

%--------------------------------------------------------------------------
% SUBFUNCTIONS
%--------------------------------------------------------------------------

function reta = R(A,x,eta)
    reta.U = qf(x.U+eta.U);
    reta.V = qf(x.V+eta.V);
    reta.W = qf(x.W+eta.W);
    reta.AU = tmprod(A,reta.U',1);
    reta.AV = tmprod(A,reta.V',2);
    reta.AW = tmprod(A,reta.W',3);
    reta.AUV = tmprod(reta.AU,reta.V',2);
    reta.AUW = tmprod(reta.AU,reta.W',3);
    reta.AVW = tmprod(reta.AV,reta.W',3);
    reta.AUVW = tmprod(reta.AUV,reta.W',3);
end

function Peta = proj(x,eta)
    Peta.U = (eye(size(x.U,1))-x.U*x.U')*eta.U;
    Peta.V = (eye(size(x.V,1))-x.V*x.V')*eta.V;
    Peta.W = (eye(size(x.W,1))-x.W*x.W')*eta.W;
end

function fval = f(A,x)
    temp = tmprod(A,{x.U,x.V,x.W},[1 2 3],'H');
    fval = -temp(:)'*temp(:);
end

function gradx = grad(x)
    [I1,R1] = size(x.U);
    [I2,R2] = size(x.V);
    [I3,R3] = size(x.W);
    gradx.U = -2*reshape(x.AVW, [I1,R2*R3])*...
                 reshape(x.AUVW,[R1,R2*R3])';
    gradx.V = -2*reshape(permute(x.AUW, [2 3 1]),[I2,R1*R3])*...
                 reshape(permute(x.AUVW,[2 3 1]),[R2,R1*R3])';
    gradx.W = -2*reshape(permute(x.AUV, [3 1 2]),[I3,R1*R2])*...
                 reshape(permute(x.AUVW,[3 1 2]),[R3,R1*R2])';
    gradx = proj(x,gradx);
end

function Heta = H(x,eta)
    [I1,R1] = size(x.U);
    [I2,R2] = size(x.V);
    [I3,R3] = size(x.W);
    Heta.U = -2*(...
        + reshape(tmprod(x.AW,eta.V',2)+...
                  tmprod(x.AV,eta.W',3),[I1,R2*R3])*...
          reshape(x.AUVW,[R1,R2*R3])'...
        + reshape(x.AVW, [I1,R2*R3])*...
          reshape(tmprod(x.AUW,eta.V',2)+...
                  tmprod(x.AUV,eta.W',3)+...
                  tmprod(x.AVW,eta.U',1),[R1,R2*R3])'...
        - eta.U*x.U'*reshape(x.AVW, [I1,R2*R3])*...
                     reshape(x.AUVW,[R1,R2*R3])');
    Heta.V = -2*(...
        + reshape(permute(tmprod(x.AU,eta.W',3)+...
                          tmprod(x.AW,eta.U',1),[2 3 1]),[I2,R1*R3])*...
          reshape(permute(x.AUVW,[2 3 1]),[R2,R1*R3])'...
        + reshape(permute(x.AUW, [2 3 1]),[I2,R1*R3])*...
          reshape(permute(tmprod(x.AUV,eta.W',3)+...
                          tmprod(x.AVW,eta.U',1)+...
                          tmprod(x.AUW,eta.V',2),[2 3 1]),[R2,R1*R3])'...
        - eta.V*x.V'*reshape(permute(x.AUW, [2 3 1]),[I2,R1*R3])*...
                     reshape(permute(x.AUVW,[2 3 1]),[R2,R1*R3])');
    Heta.W = -2*(...
        + reshape(permute(tmprod(x.AV,eta.U',1)+...
                          tmprod(x.AU,eta.V',2),[3 1 2]),[I3,R1*R2])*...
          reshape(permute(x.AUVW,[3 1 2]),[R3,R1*R2])'...
        + reshape(permute(x.AUV, [3 1 2]),[I3,R1*R2])*...
          reshape(permute(tmprod(x.AVW,eta.U',1)+...
                          tmprod(x.AUW,eta.V',2)+...
                          tmprod(x.AUV,eta.W',3),[3 1 2]),[R3,R1*R2])'...
        - eta.W*x.W'*reshape(permute(x.AUV, [3 1 2]),[I3,R1*R2])*...
                     reshape(permute(x.AUVW,[3 1 2]),[R3,R1*R2])');
    Heta = proj(x,Heta);
end

function X = qf(B)
    k = size(B,2);
    [X,R] = qr(B,0);
    diagR = diag(R);
    diagR(diagR==0) = 1;
    X = X*spdiags(sign(diagR),0,k,k);
    if rank(R) < size(B,2),
        error('qf error');
    end
end

function ez = g(~,eta,zeta)
    ez = sum(dot(eta.U,zeta.U))+ ...
         sum(dot(eta.V,zeta.V))+ ...
         sum(dot(eta.W,zeta.W));
end

function r = randT(x)  
    eta.U = randn(size(x.U));
    eta.V = randn(size(x.V));
    eta.W = randn(size(x.W));
    eta.U = eta.U / sqrt(g(x,eta,eta));
    eta.V = eta.V / sqrt(g(x,eta,eta));
    eta.W = eta.W / sqrt(g(x,eta,eta));
    r = proj(x,eta);
end

%--------------------------------------------------------------------------
% RIEMANNIAN TRUST REGION SOLVER
%--------------------------------------------------------------------------

function varargout = rtr_Gr3(fns,params)
% RTR   Riemannian Trust-Region (with tCG inner solve)
%
% This performs a retraction-based, trust-region minimization of an
% objective function on a Riemannian manifold.
%
% x = rtr(fns,params) returns the minimizer (a point on the manifold)
% [x,stats] = rtr_Gr3(fns,params) also returns some statistics from the
% algorithm
%
% The struct fns contains the function handles necessary to perform the
% optimization:
%  fns.R(x,eta) : retract tangent vector eta at T_x M to the manifold
%  fns.g(x,eta,zeta) : Riemannian metric at T_x M 
%  fns.proj(x,eta) : Project eta from nearby T_x M to T_x M (used in
%  iteration to combat round-off)
%                    Also used in debugging to check tangentiality
%                    Can be set to the identity if this is not a concern.
%  fns.f(x) : Compute the objective function at x
%  fns.fgrad(x) : Compute the gradient of f at x, returning a tangent
%  vector
%  fns.fhess(x,eta) : Apply the Hessian of f, hess_x: T_x M -> T_x M
% Optionally, three other methods may be specified:
%  fns.prec(x,eta) : apply an s.p.d. preconditioner for the inner
%  iteration, to eta, a
%                    tangent vector in T_x M
%  fns.dist(x,y) : returns the distance on the manifold between points x
%  and y
%  fns.randT(x) : returns a very small random tangent vector in T_x M
% More info:
%  fns.dist is used to measure the distance from solution and is only
%  referenced if params.xsol was specified.
%  fns.randT is used to initialize the trust-region subproblem with a
%  random vector, to encourage escape from an
%      exact critical point.
%  fns.prec is a preconditioner for the model minimization, a
%  positive-definite (undr fns.g)
%      mapping from T_x M to T_x M
%
% Parameters to the solver include:
% Required parameters:
%   params.x0        - An initial iterate. Because RTR does not know what
%   the manifold is, this is mandatory.
%   params.Delta_bar - Maximum trust-region radius. Because RTR does not
%   know anything about the manifold, this is mandatory.
%   params.Delta0    - Initial trust-region radius. Because RTR does not
%   know anything about the manifold, this is mandatory.
%   params.xsol      - For testing distance from solution using fns.dist
%   params.verbosity - Level of verbosity: 0 is silent, 1 has one-line
%   info, 2 has detailed info
%   params.debug     - Debugging tests (0 is silent, 1 has some, 2 has a
%   lot, 3 is more than you want to know)
%   params.min_outer - Minimum number of outer iterations (default: 0);
%   used only with randomization
%   params.max_outer - Maximum number of outer iterations (default: 100)
%   params.min_inner - Minimum number of inner iterations (default: 0).
%   Only in effect if using randomized initial tangent vectors.
%   params.max_inner - Maximum number of inner iterations (default: inf).
%   Recommended: dimension of manifold.
%   params.epsilon   - Outer Convergence tolerance (absolute)
%   params.kappa     - Inner kappa convergence tolerance
%   params.theta     - Inner theta convergence tolerance
%   params.rho_prime - Accept/reject ratio
%   params.testgh    - perform some simple numerical testing of gradient
%   and Hessian
%   params.useRand   - initial the trust-region solve with a random tangent
%   vector
%
% The stats output is a struct array, with fields:
%  k  - the outer iteration number for the stats
%  ng - the norm of the gradient, sqrt(g(x,gradfval,gradfval))
%  fval - the current value under the objective function
%  rho - the performance ratio for the iterate
%  time - the wallclock time for the outer iteration
%  accepted - whether the proposed iterate was accepted or not
%  numinner - the number of inner iterations used to compute the next
%  iterate
%  Delta - the trust-region radius at the outer iteration
%  dist     - the distance from the solution
%  cauchy   - whether the cauchy point was used in place of an updated
%  computed from a random initial tangent vector
%
% See also rtr, irtr

% About: RTR - Riemannian Trust-Region
% (C) 2004-2007, P.-A. Absil, C. G. Baker, K. A. Gallivan
% Florida State University
% School of Computational Science

% Modification history
% Current version: rtr_Gr3: modified by M. Ishteva, 20 Dec 2010
%   On Gr x Gr x Gr; grad and Hess have 3 components each
% Version 0.221 - CGB
%   Using preconditioner-based trust-region radius
%   Disable randomization when using preconditioner
% Version 0.2 - CGB
%   Added more commenting
%   Added preconditioning for inner iteration (fns.prec)
% Version 0.1 - CGB - Sat Nov 04 2006
%   Started from RTR03(Baker), from tCG07(Absil)

tcg_stop_reason = {'negative curvature',...
                   'exceeded trust region',...
                   'reached target residual-kappa',...
                   'reached target residual-theta',...
                   'dimension exceeded'};

if nargin < 2,
   error('Invalid arguments to rtr. See ''help rtr''');
end

% check functions handles
if   ~(isfield(fns,'R')     && isa(fcnchk(fns.R)    ,'function_handle'))...
  || ~(isfield(fns,'g')     && isa(fcnchk(fns.g)    ,'function_handle'))...
  || ~(isfield(fns,'proj')  && isa(fcnchk(fns.proj) ,'function_handle'))...
  || ~(isfield(fns,'f')     && isa(fcnchk(fns.f)    ,'function_handle'))...
  || ~(isfield(fns,'fgrad') && isa(fcnchk(fns.fgrad),'function_handle'))...
  || ~(isfield(fns,'fhess') && isa(fcnchk(fns.fhess),'function_handle')),
  error('Invalid arguments to rtr. See ''help rtr''');
end
fns.havedist = 0;
fns.haverand = 0;
fns.haveprec = 0;
if isfield(fns,'dist') && isa(fcnchk(fns.dist),'function_handle')
    fns.havedist = 1;
end
if isfield(fns,'randT') && isa(fcnchk(fns.randT),'function_handle')
    fns.haverand = 1;
else
    fns.haverand = 0;
end
if isfield(fns,'prec') && isa(fcnchk(fns.prec),'function_handle'),
    fns.haveprec = 1;
else
    fns.haveprec = 0;
end

% set default parameters
verbosity =   1;
debug     =   0;
min_inner =   0;
max_inner = inf;
min_outer =   0;
max_outer = 100;
epsilon   =   1e-6;
kappa     =   0.1;
theta     =   1.0;
rho_prime =   0.1;
%kappa_easy = .001;
testgh    = 0;
useRand  = 0;
%Delta_bar =  user must specify
%Delta0    =  user must specify
%x0        =  user must specify

% get algorithm parameters
[verbosity,params] = get_int(params,'verbosity', ...
    'Level of verbosity'                ,verbosity);
[debug    ,params] = get_int(params,'debug'    , ...
    'Debugging tests'                   ,debug);
[min_outer,params] = get_int(params,'min_outer', ...
    'Minimum number of outer iterations',min_outer,0);
[max_outer,params] = get_int(params,'max_outer', ...
    'Maximum number of outer iterations',max_outer,0);
[min_inner,params] = get_int(params,'min_inner', ...
    'Minimum number of inner iterations',min_inner,0);
[max_inner,params] = get_int(params,'max_inner', ...
    'Maximum number of inner iterations',max_inner,0);
[testgh   ,params] = get_int(params,'testgh'   , ...
    'Test gradient and Hessian'         ,testgh);
[useRand  ,params]  = get_int(params,'useRand' , ...
    'Random TR subproblem init'        ,useRand);
[x0       ,params] = get_iterate(params,'x0'   , ...
    'Initial iterate');
[xsol     ,params] = get_iterate(params,'xsol' , ...
    'Solution for testing',[]);
if isempty(xsol),
    % must have fns.dist and xsol to make fns.dist() comparisons
    fns.havedist = 0;
end
[epsilon  ,params] = get_float(params,'epsilon'  , ...
    'Outer Convergence tolerance'      ,epsilon,0);
[rho_prime,params] = get_float(params,'rho_prime', ...
    'Accept/reject parameter'          ,rho_prime,0,1);
[Delta_bar,params] = get_float(params,'Delta_bar', ...
    'Maximum trust-region radius');
[Delta0   ,params] = get_float(params,'Delta0'   , ...
    'Initial trust-region radius');
[kappa    ,params] = get_float(params,'kappa'    , ...
    'Inner kappa convergence tolerance',kappa,0,1);
[theta]            = get_float(params,'theta'    , ...
    'Inner theta convergence tolerance',theta,0);
%[kappa_easy,params] = get_float(params,'kappa_easy', ...
%    'Inner convergence tolerance'    ,kappa_easy,0,1);

if testgh && fns.haverand,
   TestGH(x0,fns);
   varargout{1} = x0;
   if nargout > 1,
      varargout{2} = struct([]);
   end
   return;
end

if fns.haveprec && useRand && fns.haverand,
   fprintf(['rtr/tCG: Cannot perform preconditioned tCG with random ' ...
            'initial vector.\n']);
   fprintf('rtr/tCG: Disabling randomization.\n');
   useRand = 0;
end

if ~fns.haverand,
   useRand = 0;
end

% ***** Initializations *****

% initialize counters/sentinals
% allocate storage for dist, counters
k          = 0;  % counter for outer (TR) iteration.
stop_outer = 0;  % stopping criterion for TR.
total_time = 0;     % total time spent in outer loop


% initialize solution and companion measures: 
% fgrad(x)
% f(x)
tic
x = x0;
fval = fns.f(x);
fgradx = fns.fgrad(x);
norm_grad = sqrt(fns.g(x,fgradx,fgradx));
this_time = toc; % not included in the total_time 

% initialize trust-region radius
Delta = Delta0;

% ** Record data:
curstat.k        = k;
curstat.ng       = norm_grad;
curstat.normS    = sqrt(-fval);
curstat.sangle   = inf;
curstat.rho      = inf;
curstat.rhonum      = 0;
curstat.rhoden      = 0;
curstat.time = this_time;
curstat.accepted = true;
curstat.numinner = 0;
curstat.Delta    = Delta;
if fns.havedist,
    curstat.dist = fns.dist(x,xsol);
end
if useRand,
    curstat.cauchy = false;
end
stats(1) = curstat;

% ** Display:
if (verbosity == 1)
   fprintf(['%3s %3s      %5s                %5s     ',...
            'f: %e   |grad|: %e\n'],...
           '   ','   ','     ','     ',fval,norm_grad);
elseif (verbosity > 1)
   fprintf('**********************************************************\n');
   fprintf('%3s %3s    k: %5s     num_inner: %5s     %s\n',...
           '','','______','______','');
   fprintf('       f(x) : %e       |grad| : %e\n',fval,norm_grad);
   fprintf('      Delta : %f\n',Delta);
   if fns.havedist,
      fprintf('       dist : %e\n',stats(end).dist);
   end
   fprintf('       Time : %f\n',this_time);
end


% **********************
% ** Start of TR loop **
% **********************
while stop_outer==0,

   % start clock for this outer iteration
   tic

   % update counter
   k = k+1;

   if (verbosity > 1) || (debug > 0),
      fprintf('*******************************************************\n');
   end

   % *************************
   % ** Begin TR Subproblem **
   % *************************
  
   % determine eta0
   if useRand,
      % random vector in T_x M
      eta = fns.randT(x);
      % must be inside trust-region
      while fns.g(x,eta,eta) > Delta^2,
         eta.U = eta.U * sqrt(sqrt(eps));
         eta.V = eta.V * sqrt(sqrt(eps));
         eta.W = eta.W * sqrt(sqrt(eps));
      end
   else
      % without randT, 0*fgradx is the only way that we 
      % know how to create a tangent vector
      eta.U = 0*fgradx.U;
      eta.V = 0*fgradx.V;
      eta.W = 0*fgradx.W;
   end
   o_eta.U = 0*eta.U;
   o_eta.V = 0*eta.V;
   o_eta.W = 0*eta.W;

   % solve TR subproblem
   [eta,numit,stop_inner] = tCG(fns,x,fgradx,eta,Delta,theta,kappa, ...
       min_inner,max_inner,useRand,debug);
   srstr = tcg_stop_reason{stop_inner};
   norm_eta = sqrt(fns.g(x,eta,eta));
   if debug > 0,
      testangle = fns.g(x,eta,fgradx) / (norm_eta*norm_grad);
   end

   % if using randomized approach, compare result with the
   % cauchy point
   % convergence proofs assume that we achieve at least the
   %   reduction of the cauchy point
   if useRand,
      used_cauchy = false;
      % check the curvature
      Hg = fns.fhess(x,fgradx);
      g_Hg = fns.g(x,fgradx,Hg);
      if g_Hg <= 0, 
         tau_c = 1;
      else
         tau_c = min( norm_grad^3/(Delta*g_Hg) , 1);
      end
      % and gen the cauchy point
      eta_c.U = -tau_c * Delta / norm_grad * fgradx.U;
      eta_c.V = -tau_c * Delta / norm_grad * fgradx.V;
      eta_c.W = -tau_c * Delta / norm_grad * fgradx.W;
      
      % now that we have computed the cauchy point in addition
      % to the returned eta, we might as well keep the better 
      % of them
      Heta = fns.fhess(x,eta);
      Heta_c = fns.fhess(x,eta_c);
      mdle  = fval + fns.g(x,fgradx,eta  ) + 1/2*fns.g(x,Heta  ,eta  );
      mdlec = fval + fns.g(x,fgradx,eta_c) + 1/2*fns.g(x,Heta_c,eta_c);
      if ( mdle > mdlec ) 
         eta = eta_c;
         norm_eta = sqrt(fns.g(x,eta,eta));
         used_cauchy = true;
      end 
   else
     %Heta = fns.fhess(x,eta);
   end 

   % compute the retraction of the proposal
   x_prop  = fns.R(x,eta);

   % compute function value of the proposal
   fval_prop = fns.f(x_prop);

   % do we accept the proposed solution or not?
   % compute the Hessian at the proposal
   Heta = fns.fhess(x,eta);

   % check the performance of the quadratic model
   rhonum = fval-fval_prop;
   rhoden = -fns.g(x,fgradx,eta) - 0.5*fns.g(x,Heta,eta);
   if debug > 0,
      if rhoden < 0,
         fprintf('Error! no model decrease!\n');
         keyboard;
      end
   end
   if debug > 0,
      fprintf('DBG:     rhonum : %e\n',rhonum);
      fprintf('DBG:     rhoden : %e\n',rhoden);
   end
   rho =   rhonum  / rhoden;
   if debug > 1,
      m = @(x,eta) fns.f(x) + fns.g(x,fns.fgrad(x),eta) + .5*fns.g(x, ...
          fns.fhess(x,eta),eta);
      actrho = (fns.f(x) - fns.f(x_prop)) / (m(x,o_eta) - m(x,eta));
      fprintf('DBG:   new f(x) : %e\n',fval_prop);
      fprintf('DBG: actual rho : %e\n',actrho);
   end
   % HEURISTIC WARNING:
   % if abs(model change) is relatively zero, we are probably near a 
   % critical point. set rho to 1.
   if abs(rhonum/fval) < sqrt(eps),
      small_rhonum = rhonum;
      rho = 1;
   else 
      small_rhonum = 0;
   end

   % choose new TR radius based on performance
   trstr = '   ';
   if rho < 1/4
      trstr = 'TR-';
      %Delta = 1/4*norm_eta;
      Delta = 1/4*Delta;
   elseif rho > 3/4 && (stop_inner == 2 || stop_inner == 1),
      trstr = 'TR+';
      %Delta = min(2*norm_eta,Delta_bar);
      Delta = min(2*Delta,Delta_bar);
   end

   % choose new iterate based on performance
   %oldgradx = fgradx;
   curstat.sangle = subspace(x.U,x_prop.U);
   if rho > rho_prime,
      accept = true;
      accstr = 'acc';
      x    = x_prop;
      fval   = fval_prop;
      fgradx = fns.fgrad(x);
      norm_grad = sqrt(fns.g(x,fgradx,fgradx));
   else
      accept = false;
      accstr = 'REJ';
   end
      
   % ** Testing for Stop Criteria
   % min_outer is the minimum number of inner iterations
   % before we can exit. this gives randomization a chance to
   % escape a saddle point.
   if norm_grad/(-fval) < epsilon && (~useRand || k > min_outer), 
       % in the original version: norm_grad instead of norm_grad/(-fval)
      stop_outer = 1;
   end

   % stop clock for this outer step
   this_time = toc;
   % update total time
   total_time = total_time + this_time;

   % ** Record data:
   curstat.k        = k;
   curstat.ng       = norm_grad;
   curstat.normS    = sqrt(-fval);
   curstat.rho      = rho;
   curstat.rhonum      = rhonum;
   curstat.rhoden      = rhoden;
   curstat.time     = this_time;
   curstat.accepted = accept;
   curstat.numinner = numit;
   curstat.Delta    = Delta;
   if fns.havedist,
      curstat.dist = fns.dist(x,xsol);
   end
   if useRand,
      curstat.cauchy = used_cauchy;
   end
   stats(end+1) = curstat;

   % ** Display:
   if verbosity == 1,
      fprintf(['%3s %3s   k: %5d     num_inner: %5d     ',...
               'f: %e   |grad|: %e   %s\n'],...
              accstr,trstr,k,numit,fval,norm_grad,srstr);
   elseif verbosity > 1,
      if useRand && used_cauchy,
         fprintf('USED CAUCHY POINT\n');
      end
      fprintf('%3s %3s    k: %5d     num_inner: %5d     %s\n',...
              accstr,trstr,k,numit,srstr);
      fprintf('       f(x) : %e     |grad| : %e\n',fval,norm_grad);
      fprintf('      Delta : %f          |eta| : %e\n',Delta,norm_eta);
      if small_rhonum ~= 0,
         fprintf('VERY SMALL rho_num: %e\n',small_rhonum);
      else
         fprintf('        rho : %e\n',rho);
      end
      if fns.havedist,
         fprintf('       dist : %e\n',stats(end).dist);
      end
      fprintf('       Time : %f\n',this_time);
   end
   if debug > 0,
      fprintf('DBG: cos ang(eta,gradf): %d\n',testangle);
      if rho==0
         keyboard;
      end
   end
   
   % stop after max_outer iterations
   if k >= max_outer,
      if (verbosity > 0),
         fprintf('\n*** timed out -- k == %d***\n',k);
      end
      stop_outer = 1;
   end 

end  % of TR loop (counter: k)
if (verbosity > 1) || (debug > 0),
   fprintf('**********************************************************\n');
end
if (verbosity > 0) || (debug > 0)
    fprintf('Total time is %f\n',total_time);
end

varargout{1} = x;
if nargout > 1,
   varargout{2} = stats;
end
end

% get_int
function [ret,opts] = get_int(opts,argname,argdesc,def,lb,ub)

    if nargin < 6
        ub = inf;
    end
    if nargin < 5,
        lb = -inf;
    end
    if nargin < 4,
        mandatory = 1;
        errstr = sprintf('%s opts.%s is mandatory, an integer in [%d,%d]',...
                         argdesc,argname,lb,ub);
    else
        mandatory = 0;
        errstr = sprintf('%s opts.%s must be an integer in [%d,%d]',...
                         argdesc,argname,lb,ub);
    end

    % Process inputs and do error-checking 
    if isfield(opts,argname)
        ret = opts.(argname);
        valid = 0;
        % check that it is an int
        if isnumeric(ret),
            ret = floor(ret);
            % check size (1 by 1) and bounds
            if isequal(size(ret),[1 1]) && lb <= ret && ret <= ub,
                valid = 1;
            end
        end
        if ~valid,
            error(errstr);
        end

        % remove field from opts
        opts = rmfield(opts,argname);
    elseif mandatory == 0,
        ret = def;
    else
        error(errorstr);
    end
end

% get_float
function [ret,opts] = get_float(opts,argname,argdesc,def,lb,ub)

    if nargin < 6
        ub = inf;
    end
    if nargin < 5,
        lb = -inf;
    end
    if nargin < 4,
        mandatory = 1;
        errstr = sprintf( ...
            '%s opts.%s is mandatory, a scalar in [%d,%d]',...
            argdesc,argname,lb,ub);
    else 
        mandatory = 0;
        errstr = sprintf( ...
            '%s opts.%s must be a scalar in [%d,%d]',...
            argdesc,argname,lb,ub);
    end

    % Process inputs and do error-checking 
    if isfield(opts,argname)
        ret = opts.(argname);
        valid = 0;
        % check that it is an int
        if isnumeric(ret),
            ret = double(ret);
            if lb <= ret && ret <= ub,
               valid = 1;
            end
        end
        if ~valid,
            error(errstr);
        end

        % remove field from opts
        opts = rmfield(opts,argname);
    elseif mandatory == 0,
        ret = def;
    else
        error(errstr);
    end
end

% get_iterate
function [ret,opts] = get_iterate(opts,argname,argdesc,def)
    if nargin < 4,
        mandatory = 1;
    else 
        mandatory = 0;
    end
    % Process inputs and do error-checking 
    if isfield(opts,argname)
        ret = opts.(argname);
        % remove field from opts
        opts = rmfield(opts,argname);
    elseif mandatory==0,
        ret = def;
    else
        error('%s opts.%s is mandatory',argdesc,argname);
    end
end

% truncated CG
function [eta,inner_it,stop_tCG] = tCG(fns,x,grad,eta,Delta,theta, ...
    kappa,min_inner,max_inner,useRand,debug)
% tCG - Truncated (Steihaug-Toint) Conjugate-Gradient
% minimize <eta,grad> + .5*<eta,Hess(eta)>
% subject to <eta,eta> <= Delta^2

   % all terms involving the trust-region radius will utilize an inner
   % product
   % w.r.t. the preconditioner; this is because the iterates grow in
   % length w.r.t. the preconditioner, guaranteeing that we will not 
   % re-enter the trust-region
   % 
   % the following recurrences for Prec-based norms and inner 
   % products come from CGT2000, pg. 205, first edition
   % below, P is the preconditioner
   % 
   % <eta_k,P*delta_k> = beta_k-1 * ( <eta_k-1,P*delta_k-1> + alpha_k-1
   % |delta_k-1|^2_P )
   % |delta_k|^2_P = <r_k,z_k> + beta_k-1^2 |delta_k-1|^2_P
   % 
   % therefore, we need to keep track of 
   % 1)   |delta_k|^2_P 
   % 2)   <eta_k,P*delta_k> = <eta_k,delta_k>_P
   % 3)   |eta_k  |^2_P
   % 
   % initial values are given by:
   %    |delta_0|_P = <r,z>
   %    |eta_0|_P   = 0
   %    <eta_0,delta_0>_P = 0
   % because we take eta_0 = 0

   if useRand, % and therefore, fns.haveprec == 0
      % eta (presumably) ~= 0 was provided by the caller
      Heta = fns.fhess(x,eta);
      r.U = grad.U+Heta.U;
      r.V = grad.V+Heta.V;
      r.W = grad.W+Heta.W;
      e_Pe = fns.g(x,eta,eta);
   else % and therefore, eta == 0
      % eta = 0*grad;
      r = grad;
      e_Pe = 0;
   end
   r_r = fns.g(x,r,r);
   norm_r = sqrt(r_r);
   norm_r0 = norm_r;

   % precondition the residual
   if fns.haveprec == 0,
      z = r;
   else
      z = fns.prec(x,r);
   end
   % compute z'*r
   z_r = fns.g(x,z,r);
   d_Pd = z_r;

   % initial search direction
   delta.U  = -z.U;
   delta.V  = -z.V;
   delta.W  = -z.W;
   if useRand, % and therefore, fns.haveprec == 0
      e_Pd = fns.g(x,eta,delta);
   else % and therefore, eta == 0
      e_Pd = 0;
   end

   % pre-assume termination b/c j == end
   stop_tCG = 5;

   % begin inner/tCG loop
   j = 0;
   for j = 1:max_inner,

      Hxd = fns.fhess(x,delta);

      % compute curvature
      d_Hd = fns.g(x,delta,Hxd);

      % DEBUGGING: check that <d,Hd> = <Hd,d>
      if debug > 1,
         Hd_d = fns.g(x,Hxd,delta);
         fprintf('DBG: |d_Hd - Hd_d| (abs/rel): %e/%e\n',abs(d_Hd-Hd_d),...
             abs((d_Hd-Hd_d)/d_Hd));
      end

      alpha = z_r/d_Hd;
      % <neweta,neweta>_P = <eta,eta>_P + 2*alpha*<eta,delta>_P +
      % alpha*alpha*<delta,delta>_P
      e_Pe_new = e_Pe + 2.0*alpha*e_Pd + alpha*alpha*d_Pd;

      if debug > 2,
         fprintf('DBG:   (r,r)  : %e\n',r_r);
         fprintf('DBG:   (d,Hd) : %e\n',d_Hd);
         fprintf('DBG:   alpha  : %e\n',alpha);
      end

      % check curvature and trust-region radius
      if d_Hd <= 0 || e_Pe_new >= Delta^2,
         % want
         %  ee = <eta,eta>_prec,x
         %  ed = <eta,delta>_prec,x
         %  dd = <delta,delta>_prec,x
         tau = (-e_Pd + sqrt(e_Pd*e_Pd + d_Pd*(Delta^2-e_Pe))) / d_Pd;
         if debug > 2,
            fprintf('DBG:     tau  : %e\n',tau);
         end
         eta.U = eta.U + tau*delta.U;
         eta.V = eta.V + tau*delta.V;
         eta.W = eta.W + tau*delta.W;
         if d_Hd <= 0,
            stop_tCG = 1;     % negative curvature
         else
            stop_tCG = 2;     % exceeded trust region
         end
         break;
      end

      % no negative curvature and eta_prop inside TR: accept it
      e_Pe = e_Pe_new;
      eta.U = eta.U + alpha*delta.U;
      eta.V = eta.V + alpha*delta.V;
      eta.W = eta.W + alpha*delta.W;

      % update the residual
      r.U = r.U + alpha*Hxd.U;
      r.V = r.V + alpha*Hxd.V;
      r.W = r.W + alpha*Hxd.W;
      % re-tangentalize r
      r = fns.proj(x,r);

      % compute new norm of r
      r_r = fns.g(x,r,r);
      norm_r = sqrt(r_r);

      % check kappa/theta stopping criterion
      if j >= min_inner && norm_r <= norm_r0*min(norm_r0^theta,kappa)
         % residual is small enough to quit
         if kappa < norm_r0^theta,
             stop_tCG = 3;  % linear convergence
         else
             stop_tCG = 4;  % superlinear convergence
         end
         break;
      end

      % precondition the residual
      if fns.haveprec == 0,
         z = r;
      else
         z = fns.prec(x,r);
      end

      % save the old z'*r
      zold_rold = z_r;
      % compute new z'*r
      z_r = fns.g(x,z,r);

      % compute new search direction
      beta = z_r/zold_rold;
      delta.U = -z.U + beta*delta.U;
      delta.V = -z.V + beta*delta.V;
      delta.W = -z.W + beta*delta.W;
      % update new P-norms and P-dots
      e_Pd = beta*(e_Pd + alpha*d_Pd);
      d_Pd = z_r + beta*beta*d_Pd;

   end  % of tCG loop
   inner_it = j;
   
end
