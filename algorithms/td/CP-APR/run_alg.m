% TD | CP-APR | PARAFAC/CP decomposition solved by Alternating Poisson Regression
% process_video('TD', 'CP-APR', 'dataset/demo.avi', 'output/demo_CP-APR.avi');

r = 10;
A = double(T);
L = double(cp_apr(T,r));
S = A - L;
