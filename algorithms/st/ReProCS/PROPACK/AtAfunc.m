function y=AtAfunc(x)
% y=AtAfunc(x)
% Testfunction  defining a linear operator applied to x. 
% Used for testing laneig.
%
%       y = A'*(A*x)

% Rasmus Munk Larsen, DAIMI, 1998


global A MxV
y = A'*(A*x);
MxV = MxV + 2;
