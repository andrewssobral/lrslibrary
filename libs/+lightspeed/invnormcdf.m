function x = invnormcdf(p)
%INVNORMCDF(P)  Normal quantile function
% X = INVNORMCDF(P) returns the P-th quantile of the standard normal distribution.
% In other words, it returns X such that P = NORMCDF(X).

x = erfinv(2*p-1)*sqrt(2);
