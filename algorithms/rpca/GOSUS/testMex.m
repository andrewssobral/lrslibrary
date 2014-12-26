
%% test shrinkageMex
clc
clear
close all

m = 200;
n = 200;
R = randn(m,n);

threshold = 0.002;
 
R(R<-threshold) = 0;
R(R>threshold) = 0;
R(abs(R) < 1e-6) = 0;

R = sparse( R );
nonZero = numel(find (R>0) )
RnonZeroRatio = nonZero/(m*n)

if RnonZeroRatio >1e-3
    return
end

display('shrinkage by Matlab')

tic
Z=R;
for i=1:size(R,2)
    r = R(:,i);
    rNorm = norm(r);

    if rNorm > 1,
        Z(:,i) = (1-1/rNorm)*r;
    else
        Z(:,i) = 0;
    end
end
toc


  display('shrinkage by mex')
  
  tic
  Z1 = shrinkageMex(R,[1,2*1e6]);
toc
  
  Diff = Z-Z1;

 diff = norm(diag(Diff))
return
% test Conjugate gradient  solver
clc
clear
close all

row = 200;
col = 100;
N = row * col;
memorycost = 8*N*N/(1024*1024*1024)
G = getGroup(row, col);
d = 5;
U =  orth(randn(N,d));
w = randn(d,1)*30;
y = abs(randn(N,1)*10);
I = eye(N);
cgParam = [0.01 10000 1e-6];

lambda =10;
rho =1.8;
Uw = U*w;
A =  G'*G +   rho*I;
tic
Ainv = inv(A);
toc
return
B  =   lambda* Uw *Uw';

% test dropping rank one matrix

x1 = A \ y;

x2 = (A+B) \ y ;

norm(x1-x2)

x1(1:10)

x2(1:10)
return



meanA = mean(abs(A(:)))
% medianA = median(abs(A(:)))

nonZeroRation=numel(find(abs(A)>5*meanA))/(N*N)
 
cgParam(1) = 5*meanA;


display('time for GPU-CG');
tic
x1 = cgCuspMex(A, y, cgParam);
toc
 
diffGPU = norm(A*x1-y)



display('time for LU');

tic
[AL, AU] = lu(A);
xtmp = AL \ y;
xtmp = AU \ xtmp;
toc
 diffLU = norm(A*xtmp-y)

 
 display('time for trucated sparse matrix')
 A1=A;
 A1(abs(A)<0.01) = 0;
%A1=sparse(A1);
tic
[AL, AU] = lu(A1);
xtmp1 = AL \ y;
xtmp1 = AU \ xtmp1;
toc
 diffLUtrucated = norm(A*xtmp1-y)

return

%% Warning: The device selected (device 1, "GeForce GT
% 330") does not have sufficient compute capability to be
% used. Compute capability 1.3 (or greater) is required,
% the selected device has compute capability 1.2. 
display ('time for GPU');
tic
A_gpu = gpuArray(A);
y_gpu = gpuArray(y);
x_gpu = A \ y;
x= gather(x_gpu);
toc
return

load ../test/testGroupSparisty/testData/linearSystem.mat
meanA = mean(abs(A(:)))
% medianA = median(abs(A(:)))
cgParam(1) = 5*meanA;
y = b;
display('time for GPU-CG');
tic
x1 = cgCuspMex(A, y, cgParam);
toc



