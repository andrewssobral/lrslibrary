%mxn data matrix M with rank r
m = 200;
n = 300;
r = 3;

U0 = rand(m,r);
V0 = rand(r,n);

W = double(rand(m,n) > 0.3);

%missing data ratio
missing_ratio = 1 - sum(sum(W))/(m*n);

%ground truth 
M0 = U0*V0;

%adding Gaussian noise
M = M0 + rand(m,n)*0.01;

%outliers
G = double(rand(m,n) > 0.9);

%outlier ratio
outlier_ratio = sum(sum(G))/(m*n);

%adding outliers with fixed magnatude of 1
M = M + G*1;

%call RegL1-ALM, 
%NOTE: for randomly generated data, it is usually enough to set 
%rho = 1.2, maxIterIN = 1, 
[M_est U_est V_est L1_error]=RobustApproximation_M_UV_TraceNormReg(M,W,r,1e-3,1.2,1,0);

%%
W = ones(size(M));
r = 1;
[M_est,U_est,V_est,L1_error] = ...
  RobustApproximation_M_UV_TraceNormReg(M,W,r,1e-3,1.2,1,0);

%%
L = M_est;
S = M - L;
show_2dvideo(S,video.height,video.width);

%%

%The average L1-norm error over observed entries w.r.t the noisy M
Ave_L1_error_M = sum(sum(abs(W.*(M-M_est))))/(sum(sum(W)));

%The average L1-norm error over observed entries w.r.t the ground-truth M0
Ave_L1_error_M0 = sum(sum(abs((M0-M_est))))/(m*n);

disp('Missing data ratio:');
fprintf('%f%%\n', missing_ratio*100);
    
disp('Outlier ratio');
fprintf('%f%%\n', outlier_ratio*100);
    
if m > 5 && n > 5    
    disp('Ground Truth M0(1:5,1:5):');
    disp(M0(1:5,1:5));
    
    disp('Noisy M(1:5,1:5):');
    disp(M(1:5,1:5));
    
    disp('Estimated M_est(1:5,1:5):');
    disp(M_est(1:5,1:5));
end

%draw error w.r.t M
disp('A figure with two straigh lines at 1 and 0, respectively, means GOOD: ');
Error = abs(M - M_est);
plot(Error(:),'b.');

    

