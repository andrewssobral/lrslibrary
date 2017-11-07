function Y = ifWHtrans(X)
% ifWHtrans computes fast inverse discrete Walsh-Hadamard transform 
% with sequency order.
%
% Since the forward and inverse transforms are exactly identical
% operations, fastWHtrans.mexw32 is used to perform inverse transform.
%
% written by: Chengbo Li
% Computational and Applied Mathematics Department, Rice University.
% 02/15/2010

Y = fWHtrans(X);
% Perform scaling
[m, n] = size(Y);
Y = Y .* m;
