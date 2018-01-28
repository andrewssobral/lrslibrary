function y = Cfunc(x)
% y=Cfunc(x)
% Testfunction  defining a linear operator applied to x. 
% Used for testing laneig.
%
%       y =     [ 0  A ]  * x
%               [ A' 0 ]

% Rasmus Munk Larsen, DAIMI, 1998


global A MxV
[m n] = size(A);
y = [A*x(m+1:end,:); A'*x(1:m,:)];
MxV = MxV + 2;