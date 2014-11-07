%              ==== VERY IMPORTANT ====
% This code needs the support of Tensor Toolbox developed by Tamara Kolda 
%  which is available at:
%         http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html
%
clear;
clc;
I=[50,50,50];
R=[5,6,7];
N=numel(I);

% Generate data;
A=cell(N,1);
for n=1:N
    A{n}=rand(I(n),R(n));
end
Y=ttensor(tensor(rand(R)),A);
Y=tensor(Y);


opts=struct('NumOfComp',R,'nlssolver','hals','maxiter',100,'maxiniter',20,'tdalgFile','call_tucker_als_opts.mat');
tic;
[Ydec]=lraNTD_ANLS(Y,opts);
toc;
fprintf('Complete. Fit=%f\n',fitness(Y,Ydec));
