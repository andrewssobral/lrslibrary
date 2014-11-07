function [] = demo()
A = randn(100,200);
X = randn(100,100);
lambda = 0.01;

disp('solve min |Z|_* + lambda |E|_21 s.t. X = AZ + E by exact ALM ...');

tic;
[Z1,E1] = solve_lrr(X,A,lambda,0,0);
obj1 = sum(svd(Z1)) + lambda*sum(sqrt(sum(E1.^2,1)));
toc;

disp(['objective value=' num2str(obj1)]);

disp('solve min |Z|_* + lambda |E|_21 s.t. X = AZ + E by inexact ALM ...');

tic;
[Z2,E2] = solve_lrr(X,A,lambda,0,1);
obj2 = sum(svd(Z2)) + lambda*sum(sqrt(sum(E2.^2,1)));
toc;
disp(['objective value=' num2str(obj2)]);

diff = max(max(abs(Z1 - Z2)));

warning(['difference of the solution found by those two approaches: |Z1 - Z2|_inf=' num2str(diff)]);

disp('solve min |Z|_* + lambda |E|_1 s.t. X = AZ + E by exact ALM ...');
tic;
[Z1,E1] = solve_lrr(X,A,lambda,1,0);
obj1 = sum(svd(Z1)) + lambda*sum(sqrt(sum(E1.^2,1)));
toc;

disp(['objective value=' num2str(obj1)]);

disp('solve min |Z|_* + lambda |E|_1 s.t. X = AZ + E by inexact ALM ...');
tic;
[Z2,E2] = solve_lrr(X,A,lambda,1,1);
obj2 = sum(svd(Z2)) + lambda*sum(sqrt(sum(E2.^2,1)));
toc;
disp(['objective value=' num2str(obj2)]);

diff = max(max(abs(Z1 - Z2)));

warning(['difference of the solution found by those two approaches: |Z1 - Z2|_inf=' num2str(diff) ]);
