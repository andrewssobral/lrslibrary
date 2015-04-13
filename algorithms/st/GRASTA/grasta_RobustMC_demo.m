%%
% % % % % % % % % % % % % % % % % % % % % %
% Sample code to run GRASTA: Robust Matrix Completion
% % % % % % % % % % % % % % % % % % % % % %
%
% Jun He and Laura Balzano, Jan. 2012
% 


%% First we generate the data matrix and incomplete sample vector.
clc; close all;
clear all;

grasta_path; % Add search path
%%
% s = RandStream('swb2712','Seed',1020);
% RandStream.setDefaultStream(s);

nMonteCarlo = 10;

for iTest = 1: nMonteCarlo,
    
outlierFac = 0.1;
SAMPLING   = 0.3;
noiseFac   = 1 * 1e-6;

% Number of rows and columns
numr = 1000;
numc = 1000;

probSize = [numr,numc];
% Rank of the underlying matrix.
truerank = 10;

% Size of vectorized matrix
N = numr*numc;
% Number of samples that we will reveal.
M = round(SAMPLING * N);

% The left and right factors which make up our true data matrix Y.
YL = randn(numr,truerank); 
YR = randn(numc,truerank); 

% Select a random set of M entries of Y.
p = randperm(N);
idx = p(1:M);
clear p;

[I,J] = ind2sub([numr,numc],idx(1:M));
[J, inxs]=sort(J'); I=I(inxs)';

% S denotes the values of Y at the locations indexed by I and J.
S = sum(YL(I,:).*YR(J,:),2);

% Add Gaussian noise.
noise = noiseFac*max(S)*randn(size(S));
S = S + noise;

% Add sparse outliers
outlier_magnitude = 0.01*max(abs(S));
idx = randperm(length(S));
sparseIdx = idx(1:ceil(outlierFac*length(S)))';
Outlier_part = outlier_magnitude * randn(size(sparseIdx));

S(sparseIdx) = S(sparseIdx) + Outlier_part;

% Now we set parameters for the algorithm.
% We set the number of cycles and put the necessary parameters into OPTIONS

maxCycles                   = 10;    % the max cycles of robust mc
OPTIONS.QUIET               = 1;     % suppress the debug information

OPTIONS.MAX_LEVEL           = 20;    % For multi-level step-size,
OPTIONS.MAX_MU              = 15;    % For multi-level step-size
OPTIONS.MIN_MU              = 1;     % For multi-level step-size

OPTIONS.DIM_M               = numr;  % your data's ambient dimension
OPTIONS.RANK                = truerank; % give your estimated rank

OPTIONS.ITER_MIN            = 20;    % the min iteration allowed for ADMM at the beginning
OPTIONS.ITER_MAX            = 20;    % the max iteration allowed for ADMM
OPTIONS.rho                 = 2;   % ADMM penalty parameter for acclerated convergence
OPTIONS.TOL                 = 1e-8;   % ADMM convergence tolerance

OPTIONS.USE_MEX             = 1;     % If you do not have the mex-version of Alg 2
                                     % please set Use_mex = 0.
                                     
CONVERGE_LEVLE              = 20;    % If status.level >= CONVERGE_LEVLE, robust mc converges

% % % % % % % % % % % % % % % % % % % % %
% Now run robust matrix completion.
t0 = clock;

[Usg, Vsg, Osg] = grasta_mc(I,J,S,numr,numc,maxCycles,CONVERGE_LEVLE,OPTIONS);

StochTime = etime(clock,t0);

%
% % % % % % % % % % % % % % % % % % % % % %
% Compute errors.
% % % % % % % % % % % % % % % % % % % % % %

% Select a random set of entries of Y
p = randperm(N);
idx2 = p(1:M);
clear p;

ITest = mod(idx2(1:M)-1,numr)'+1;  % row values
JTest = floor((idx2(1:M)-1)/numr)'+1; % column values

% STest denotes the values of Y at the locations indexed by ITest and JTest
STest =  sum(YL(ITest,:).*YR(JTest,:),2);


% % % % % % % % % % % % % % % % % % % % % %
% Error for Robust MC

dev = sum(Usg(ITest,:).*Vsg(JTest,:),2) - STest;
RelDevStochG = norm(dev,'fro') / norm(STest,'fro');
fprintf(1,'%d/%d, GRASTA robust mc [%d*%d, rank=%d, outlier=%.1f%%, sampling=%.1f%%]: Rel.err on test = %7.2e in %4.2f seconds\n',...
    iTest, nMonteCarlo,numr,numc, truerank, outlierFac*100, SAMPLING*100,RelDevStochG, StochTime);


end
%%
% % % % % % % % % % % % % % % % % % % % % %
% Random select one column to show the recovery result.
% % % % % % % % % % % % % % % % % % % % % %
p = randperm(numc); iSelect = p(1);

S_free = YL*YR';
Col_data = S_free(:,iSelect);

Col_max = max(abs(Col_data));
idx = randperm(numr); sparseIdx = idx(1:ceil(outlierFac*numr))';
Col_outlier = zeros(numr,1);
Col_outlier(sparseIdx) = 1 * Col_max * randn(size(sparseIdx));

Col_noise = noiseFac * Col_max * randn(numr,1);

Col_data = Col_data + Col_outlier + Col_noise;


OPTS2.MAX_ITER      = 30;
OPTS2.TOL           = 1e-8;

if OPTIONS.USE_MEX,
    [s, w, ~] = mex_srp(Usg, Col_data, OPTS2);
else
    [s, w, ~, ~] = admm_srp(Usg, Col_data, OPTS2);
end

Col_outlier_hat = s;
Col_hat = Usg * w;

Col_err = norm(Col_hat - S_free(:, iSelect))/norm(S_free(:, iSelect));
Col_outlier_err = norm(Col_outlier_hat - Col_outlier)/norm(Col_outlier);
%%%%%%%%%%%%%%
figure; 
fontsize = 16;
subplot(3,1,1);
plot(Col_outlier,'r');ylabel('$s_0$','Interpreter','latex','fontsize',fontsize); 
title('Outliers','fontsize',fontsize); axis tight;

subplot(3,1,2);
plot(Col_outlier_hat,'b');ylabel('$\hat{s}$','Interpreter','latex','fontsize',fontsize); 
title('Separated outliers','fontsize',fontsize); axis tight;

subplot(3,1,3);
plot(abs(Col_outlier-Col_outlier_hat),'b'); 
ylabel('$|\hat{s}$ - $s_0|$','Interpreter','latex','fontsize',fontsize);
xlabel('Column components','fontsize',fontsize);
title(['Oulter rel-err ' num2str(Col_outlier_err)],'fontsize',fontsize); axis tight;

%%%%%%%%%%%%%%
figure;
subplot(3,1,1);
plot(S_free(:,iSelect),'r');ylabel('$v_0$','Interpreter','latex','fontsize',fontsize); 
title('Clean column','fontsize',fontsize); axis tight;

subplot(3,1,2);
plot(Col_hat,'b');ylabel('$\hat{v}$','Interpreter','latex','fontsize',fontsize); 
title('Recovered column','fontsize',fontsize); axis tight;

subplot(3,1,3);
plot(abs(Col_hat-S_free(:,iSelect)),'b'); 
ylabel('$|\hat{v}$ - $v_0|$','Interpreter','latex','fontsize',fontsize);
xlabel('Column components','fontsize',fontsize);
title(['Column rel-err ' num2str(Col_err)],'fontsize',fontsize); axis tight;
