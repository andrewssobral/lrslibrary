% TD | t-SVD | Tensor SVD in Fourrier Domain (Zhang et al. 2013)
% process_video('TD', 't-SVD', 'dataset/demo.avi', 'output/demo_t-SVD.avi');

A = double(T);
[U,S,V] = tensor_t_svd(A);
[C] = tensor_product(U,S);
[L] = tensor_product(C,tensor_transpose(V));
S = A - L;
