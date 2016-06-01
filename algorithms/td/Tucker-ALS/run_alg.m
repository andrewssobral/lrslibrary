% TD | Tucker-ALS | Tucker Decomposition solved by Alternating Least Squares
% process_video('TD', 'Tucker-ALS', 'dataset/demo.avi', 'output/demo_Tucker-ALS.avi');

r = [size(T,1) size(T,2) 1];
A = double(T);
L = double(tucker_als(T,r));
S = A - L;
