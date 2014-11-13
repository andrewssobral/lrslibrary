%{
Recreates some tests from the conference paper
"A variational approach to stable principal component pursuit", 
A. Aravkin, S. Becker, V. Cevher, P. Olsen
UAI 2014
http://arxiv.org/abs/1406.1089


This script runs the code. There are several "fancy" things:

1) it calls loadSyntheticProblem( problemName )
     This lets you load problems that have pre-computed reference solutions
     It also lets you load problems and the corresponding parameters
       for all variants (objectives vs. constraints). It does
       this by running the Lagrangian version, from which we can
       calculate the other parameters.

2) it can run the software NSA or PSPG or ASALM for comparison, if these
    codes are visible on the Matlab path
    (Note: it runs modified versions of these, in order to collect
     errors and such, plus minor improvements, but we do not distribute
     these codes.

3) it has a fancy error function that not only records errors per
   iteration, but also records at what time the error was collected,
   and it can subtract off the time taken to collect the error,
   and it can also make a program quit if it exceeds a given time limit.
  (if we exceed the time limit, it throws the error:
    errorFunction:timedOut
   This is the reason for the "try"/"catch" statements below)

The paper used some synthetic problems taken from the NSA paper,
and since these were not in Lagrangian format, we don't collect all
parameters, so we cannot test all variants

Stephen Becker, November 12 2014
%}

tolFactor   = 1;
maxIts      = 50;
maxTime     = 20; % in seconds
doNSA              = true && exist('nsa','file');
doPSPG             = true && exist('pspg','file');
doASALM            = true && exist('asalm','file');

% Pick a problem:
% The "maxTime" is a time limit (in seconds) for each solver
problemString = 'medium_exponentialNoise_v2'; problem_type=0;tolFactor = 1e-3;maxIts = 100;maxTime=20;
% problemString = 'nsa1500_SNR45_problem1_run1'; problem_type=1;maxTime = 60; tolFactor=1e-2;% 
% problemString = 'nsa1500_SNR45_problem2_run1'; problem_type=2;maxTime = 60; tolFactor=1e-2;% 



[Y,L,S,params] = loadSyntheticProblem(problemString);
rng(100); % this affects the randomized SVDs
normL = norm(L(:)); normS = norm(S(:));
errFcn  = @(LL,SS) norm(LL-L,'fro')/normL + norm(SS-S,'fro')/normS;
erF     = @(varargin)errorFunction( errFcn, varargin{:} );

%% Solve via NSA
if doNSA
    tol     = 1e-2*tolFactor;
    t1      = tic;
    errorFunction( {t1, maxTime } ); % only allow 20 seconds
    try
        [L,S] = nsa(Y,-params.epsilon,tol,0,[],[],params.lambdaSum, erF  );
        tm      = toc(t1);
        fprintf('Via NSA, time is %.2e, error is %.2e\n', tm, errFcn(L,S) );
    catch thrownError
        if ~strcmpi(thrownError.identifier,'errorFunction:timedOut')
            rethrow(thrownError)
        else
            disp('NSA timed out, but error were still recorded');
        end
    end
    [errHist_NSA, times_NSA]     = errorFunction();
end

%% Solve via ASALM, constrained version
if doASALM 
    %   If 'minType' is 'unc', solves Lagrangian version
    Omega = (1:numel(Y))';
    opts = struct('minType','con','tol',1e-3*tolFactor,'print',1,'record_res',true);
    if problem_type == 1
        opts.beta = 0.1/mean(abs(Y(Omega)));
    else
        opts.beta = 0.15/mean(abs(Y(Omega)));
    end
    if numel(Y) <= 50^2
        opts.SVD = true; opts.PROPACK = false;
    end
    
    t1  = tic;
    errorFunction( {t1, maxTime } ); % only allow 20 seconds
    try
        out = asalm(size(Y),Omega,Y(Omega),params.lambdaSum,params.epsilon,opts,[],[],[],[],erF); %,optimal_X,optimal_S);
        tm      = toc(t1);
        L = out.LowRank; S = out.Sparse;
        fprintf('Via ASALM, constrained, time is %.2e, error is %.2e\n', tm, errFcn(L,S) );
    catch thrownError
        if ~strcmpi(thrownError.identifier,'errorFunction:timedOut')
            rethrow(thrownError)
        else
            disp('ASALM timed out, but error were still recorded');
        end
    end
    [errHist_ASALM,times_ASALM]   = errorFunction();
end
%% Solve via ASALM, unconstrained
if doASALM
    if ~isempty( params.lambdaS )
        Omega = (1:numel(Y))';
        opts = struct('minType','unc','tol',1e-3*tolFactor,'print',1,'record_res',true);
        if problem_type == 1
            opts.beta = 0.1/mean(abs(Y(Omega)));
        else
            opts.beta = 0.15/mean(abs(Y(Omega)));
        end
        if numel(Y) <= 50^2
            opts.SVD = true; opts.PROPACK = false;
        end
        
        t1  = tic;
        errorFunction( {t1, maxTime } ); % only allow 20 seconds
        try
            out = asalm(size(Y),Omega,Y(Omega),params.lambdaSum,params.epsilon,opts,[],[],[],[],erF); %,optimal_X,optimal_S);
            tm      = toc(t1);
            L = out.LowRank; S = out.Sparse;
            fprintf('Via ASALM, unconstrained, time is %.2e, error is %.2e\n', tm, errFcn(L,S) );
        catch thrownError
            if ~strcmpi(thrownError.identifier,'errorFunction:timedOut')
                rethrow(thrownError)
            else
                disp('ASALM timed out, but error were still recorded');
            end
        end
        [errHist_ASALM2,times_ASALM2]   = errorFunction();
    else
        errHist_ASALM2 = [];
    end
end
%% Solve via PSPG...
if doPSPG
  tol     = 1e-3*tolFactor;
  mu = 1.25/norm(Y); % their default
  t1  = tic;
  errorFunction( {t1, maxTime } ); 
  try
      [L,S,out]=pspg(Y,-params.epsilon,tol,[],[],params.lambdaSum,erF,mu,maxIts);
      tm      = toc(t1);
      fprintf('Via PSPG, constrained, time is %.2e, error is %.2e\n', tm, errFcn(L,S) );
  catch thrownError
      if ~strcmpi(thrownError.identifier,'errorFunction:timedOut')
          rethrow(thrownError)
      else
          disp('PSPG timed out, but error were still recorded');
      end
  end
  [errHist_PSPG,times_PSPG]   = errorFunction();
end

%% Test our dual version by itself, max version (not doing SPGL1 stuff)
opts = struct('printEvery',5,'tol',1e-4*tolFactor,'max',true,'maxIts',maxIts);
opts.errFcn     = erF;
opts.quasiNewton  = true; opts.FISTA=false; opts.BB = false;
t1  = tic;
errorFunction( {t1, maxTime } ); 
try
    [L,S] = solver_RPCA_constrained(Y,params.lambdaMax, params.tauMax, [], opts);
    tm      = toc(t1);
    fprintf('Via our "dual-max" code, constrained, time is %.2e, error is %.2e\n', tm, errFcn(L,S) );
catch thrownError
    if ~strcmpi(thrownError.identifier,'errorFunction:timedOut')
        rethrow(thrownError)
    else
        disp('our "dual-max" code timed out, but error were still recorded');
    end
end
[errHist_max,times_max]   = errorFunction();

%% Test our dual version by itself, sum version (not doing SPGL1 stuff)
opts = struct('printEvery',5,'tol',1e-5*tolFactor,'sum',true,'maxIts',maxIts);
opts.restart    = 100;
opts.errFcn     = erF;
opts.BB         = true;
opts.FISTA      = false;

t1  = tic;
errorFunction( {t1, maxTime } );
try
    [L,S] = solver_RPCA_constrained(Y,params.lambdaSum, params.tauSum, [], opts);
    tm      = toc(t1);
    fprintf('Via our "dual-sum" code, constrained, time is %.2e, error is %.2e\n', tm, errFcn(L,S) );
catch thrownError
    if ~strcmpi(thrownError.identifier,'errorFunction:timedOut')
        rethrow(thrownError)
    else
        disp('our "dual-sum" code timed out, but error were still recorded');
    end
end
[errHist_sum,times_sum]   = errorFunction();

%% Our Lagrangian version
if ~isempty( params.lambdaS )
    opts = struct('printEvery',1,'tol',1e-4*tolFactor,'maxIts',maxIts,...
        'errFcn',erF,'quasiNewton',true,'FISTA',false,'BB',false,...
        'SVDstyle',4,'SVDnPower',2);
    t1  = tic;
    errorFunction( {t1, maxTime } ); % only allow 20 seconds
    try
        [L,S] = solver_RPCA_Lagrangian(Y,params.lambdaL, params.lambdaS, [], opts);
        tm      = toc(t1);
        fprintf('Via our "Lagrangian" code, constrained, time is %.2e, error is %.2e\n', tm, errFcn(L,S) );
    catch thrownError
        if ~strcmpi(thrownError.identifier,'errorFunction:timedOut')
            rethrow(thrownError)
        else
            disp('our "Lagrangian" code timed out, but error were still recorded');
        end
    end
    [errHist_Lag,times_Lag]   = errorFunction();
else
    errHist_Lag = [];
end

%% Our Flip-Flop max version
opts = struct('printEvery',500,'tol',1e-4*tolFactor,'max',true,'maxIts',maxIts);
opts.errFcn     = erF;
opts.quasiNewton  = true; opts.FISTA=false; opts.BB = false;
opts.SPGL1_tol  = 1e-2*tolFactor;
t1  = tic;
errorFunction( {t1, maxTime } ); 
try
    [L,S] = solver_RPCA_SPGL1(Y,params.lambdaMax, params.epsilon, [], opts);
    tm      = toc(t1);
    fprintf('Via our "Flip-Flop max" code, constrained, time is %.2e, error is %.2e\n', tm, errFcn(L,S) );
catch thrownError
    if ~strcmpi(thrownError.identifier,'errorFunction:timedOut')
        rethrow(thrownError)
    else
        disp('our "Flip-Flop max" code timed out, but error were still recorded');
    end
end
[errHist_FlipMax,times_FlipMax]   = errorFunction();

%% Our Flip-Flop sum version
opts = struct('printEvery',500,'tol',1e-4*tolFactor,'sum',true,'maxIts',maxIts);
opts.errFcn     = erF;
opts.BB         = true;
opts.FISTA      = false;
opts.SVDstyle   = 4;
opts.SPGL1_tol  = 1e-2*tolFactor;
t1  = tic;
errorFunction( {t1, maxTime } );
try
    [L,S] = solver_RPCA_SPGL1(Y,params.lambdaSum, params.epsilon, [], opts);
    tm      = toc(t1);
    fprintf('Via our "Flip-Flop sum" code, constrained, time is %.2e, error is %.2e\n', tm, errFcn(L,S) );
catch thrownError
    if ~strcmpi(thrownError.identifier,'errorFunction:timedOut')
        rethrow(thrownError)
    else
        disp('our "Flip-Flop sum" code timed out, but error were still recorded');
    end
end
[errHist_FlipSum,times_FlipSum]   = errorFunction();

%% plot all of them
figure(2); clf;
holdMarker('reset');
fontsize = 14;
style = {'-','linewidth',2,'markersize',8};

% Plot error vs. time
plotF   = @semilogy;xString = 'time (s)';
% If you want to plot error vs. iterations, not time, use the following:
% plotF   = @(x,varargin) semilogy(1:length(x), varargin{:} ); xString = 'iterations';

lString = {};

plotF( times_max, errHist_max , style{:} );lString{end+1} = 'this paper* (flip-SPCP-max, QN)';
hold all
plotF( times_sum, errHist_sum , style{:} );lString{end+1} = 'this paper (flip-SPCP-sum)';
if ~isempty( errHist_Lag )
    plotF( times_Lag, errHist_Lag , style{:} );lString{end+1} = 'this paper* (lag-SPCP, QN)';
end
plotF( times_FlipMax, errHist_FlipMax, style{:} ); lString{end+1} = 'this paper* (SPCP-max, QN)';
plotF( times_FlipSum, errHist_FlipSum, style{:} ); lString{end+1} = 'this paper (SPCP-sum)';

if doNSA 
    plotF( times_NSA, errHist_NSA, style{:} ); lString{end+1} = 'NSA* (SPCP-sum)';
end
if doPSPG
    plotF( times_PSPG, errHist_PSPG , style{:} );lString{end+1} = 'PSPG* (SPCP-sum)';
end
if doASALM
    plotF( times_ASALM, errHist_ASALM , style{:} );lString{end+1} = 'ASALM* (SPCP-sum)';
    if ~isempty( errHist_ASALM2 )
        plotF( times_ASALM2, errHist_ASALM2 , style{:} );lString{end+1} = 'ASALM* (lag-SPCP)';
    end
end
% The ^*  denotes fastSVD, and the ^{QN} denotes quasi-Newton

holdMarker(0);
lh = legend(lString); set(lh,'fontsize',fontsize);
refresh
ylabel('Error','fontsize',fontsize);
xlabel(xString,'fontsize',fontsize);
switch problem_type
    case 0
        ylim([1e-4,Inf]);
    case 1
        ylim([5e-3,3]);
    case 2
        ylim([1e-3,5]);
        set(lh,'location','east');
        ps = get(gcf,'position');
        set(gcf,'position',[ps(1:2), 600,400] );
end