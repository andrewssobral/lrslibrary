%% Shifted Symmetric Higher-Order Power Method (SSHOPM)

%% Data tensor 
% From Example 1 in E. Kofidis and P. A. Regalia, On the best rank-1
% approximation of higher-order supersymmetric tensors, SIAM J. Matrix
% Anal. Appl., 23 (2002), pp. 863–884, DOI: 10.1137/S0895479801387413.
A = tenzeros([3 3 3 3]);
A(perms([1 1 1 1])) = 0.2883;
A(perms([1 1 1 2])) = -0.0031;
A(perms([1 1 1 3])) = 0.1973;
A(perms([1 1 2 2])) = -0.2485;
A(perms([1 1 2 3])) = -0.2939;
A(perms([1 1 3 3])) = 0.3847;
A(perms([1 2 2 2])) = 0.2972;
A(perms([1 2 2 3])) = 0.1862;
A(perms([1 2 3 3])) = 0.0919;
A(perms([1 3 3 3])) = -0.3619;
A(perms([2 2 2 2])) = 0.1241;
A(perms([2 2 2 3])) = -0.3420;
A(perms([2 2 3 3])) = 0.2127;
A(perms([2 3 3 3])) = 0.2727;
A(perms([3 3 3 3])) = -0.3054;

%% Call SSHOPM with no shift
% The method with no shift will fail to converge.
[lambda, x, flag, it, ~, trace] = sshopm(A, 'MaxIts', 100);
plot(0:it, trace(1:it+1),'b.-');

%% Call SSHOPM with shift

[lambda, x, flag, it, ~, trace] = sshopm(A, 'MaxIts', 100, 'Shift', 1);
plot(0:it, trace(1:it+1),'b.-');
