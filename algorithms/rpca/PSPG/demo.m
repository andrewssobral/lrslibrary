clc
clear all;
rng(13) %set the seed for random number generator
addpath ./PROPACK_SVT/
tol = 0.05;
SNR = 45; %Signal to Noise Ratio of the data matrix
n1 = 500; n2= 500;
p_ratio = 0.1;
r_ratio = 0.1;

%Data matrix generation
p = ceil(n1*n2*p_ratio);
r = ceil(min(n1,n2)*r_ratio);
mg = sqrt(8*r/pi);
pwr = r+p/(n1*n2)*mg^2/3;
stdev=sqrt(pwr/10^(SNR/10));
optimal_S = zeros(n1,n2);
optimal_X = randn(n1,r)*randn(r,n2);
noise = stdev*randn(n1,n2);
OMEGA = randperm(n1*n2);
OMEGA = OMEGA(1:p);
optimal_S(OMEGA) = -mg+2*mg*rand(1,p);
D = optimal_X+optimal_S+noise;

%Call PSPG to solve SPCP problem
[X,S,out] = pspg(D,stdev,tol,optimal_X,optimal_S);