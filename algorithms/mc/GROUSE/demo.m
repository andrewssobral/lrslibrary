clc;clear all;
%%
% % % % % % % % % % % % % % % % % % % % % %
% Sample code to run GROUSE: Stochastic Matrix Completion
% % % % % % % % % % % % % % % % % % % % % %
%
% This helper function shows how to set up a matrix completion problem and how to call
% grouse for matrix completion. If you would like code that runs to a
% particular residual error, or code that has diagnostics built in, please
% email Laura Balzano: sunbeam@ece.wisc.edu.
%
% Ben Recht and Laura Balzano, February 2010
%
%
%% First we generate the data matrix and incomplete sample vector.

% Number of rows and columns
numr = 3000;
numc = 10000;
probSize = [numr,numc];
% Rank of the underlying matrix.
truerank = 5;

% Size of vectorized matrix
N = numr*numc;
% Number of samples that we will reveal.
M = floor(.026*N);

% The left and right factors which make up our true data matrix Y.
YL = randn(numr,truerank);
YR = randn(numc,truerank);

% Select a random set of M entries of Y.
idx = randperm(N);
idx = idx(1:M);

[I,J] = ind2sub([numr,numc],idx);
[J, inxs]=sort(J'); I=I(inxs)';

% Values of Y at the locations indexed by I and J.
S = sum(YL(I,:).*YR(J,:),2);
S_noiseFree = S;

% Add noise.  
noiseFac = 0;
noise = noiseFac*randn(size(S));
S = S + noise;


%% Now we set parameters for the algorithm.

% We set an upper bound for the rank of the matrix and the number of cycles
% or passes over the data that we want to execute. We also set the gradient
% step size. There is an order of magnitude of step sizes for which the
% algorithm converges to the nosie level.
maxrank = truerank;
maxCycles = 3;
step_size = 0.1;

%% % % % % % % % % % % % % % % % % % % % %
% Now run GROUSE.
t0 = clock;

[Usg, Vsg, err_reg] = grouse(I,J,S,numr,numc,maxrank,step_size,maxCycles);

StochTime = etime(clock,t0);

%%
% % % % % % % % % % % % % % % % % % % % % %
% Compute errors.  Assuming that Y = YL*YR' fits in memory
% % % % % % % % % % % % % % % % % % % % % %

fprintf(1,'\n Results: \n');

% % % % % % % % % % % % % % % % % % % % % %
% Error for Stochastic gradient

dev = norm(Usg*Vsg' - YL*YR','fro');
RelDevStochG = norm(dev,'fro') / norm(YL*YR','fro');
fprintf(1,'stoch grad rel.err = %7.2e in %4.2f seconds\n', RelDevStochG, StochTime);


