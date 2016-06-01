% TD | CP-ALS | PARAFAC/CP decomposition solved by Alternating Least Squares
% process_video('TD', 'CP-ALS', 'dataset/demo.avi', 'output/demo_CP-ALS.avi');

r = 10;
A = double(T);
L = double(cp_als(T,r,'dimorder',[3 2 1]));
S = A - L;
