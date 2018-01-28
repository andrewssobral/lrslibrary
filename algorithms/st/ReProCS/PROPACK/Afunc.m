function y=Afunc(x)
% y=Afunc(x)
% Testfunction returning a linear operator applied to x.
% Used for testing lansvd.
%
% y = A'*x

% Rasmus Munk Larsen, DAIMI, 1998

global A MxV
y = A*x;
MxV = MxV + 1;