function [] = test()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
randn('state',1212412414424234324);
rand('state',1212412414424234324);
dim_ambient = 1000;
dim_r = 10;
M = randn(dim_ambient,dim_r);
N = randn(dim_r,dim_ambient);
D0 = M*N/dim_r;

E0 = sign(randn(dim_ambient,dim_ambient));
inds = rand(dim_ambient)<0.7;
E0(inds) = 0;

D = D0 + E0;

D_hat = as_rpca(D,0.05,5,1.3);

error = max(max(abs(D0 - D_hat)))./max(max(abs(D0)));
disp(['recover error=' num2str(error)]);

D_hat = as_rpca(D,0.05,10,1.3);

error = max(max(abs(D0 - D_hat)))./max(max(abs(D0)));
disp(['recover error=' num2str(error)]);

D_hat = as_rpca(D,0.05,15,1.3);

error = max(max(abs(D0 - D_hat)))./max(max(abs(D0)));
disp(['recover error=' num2str(error)]);
end

