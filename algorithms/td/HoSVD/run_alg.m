% TD | HoSVD | High-order singular value decomposition (Tucker decomposition)
% process_video('TD', 'HoSVD', 'dataset/demo.avi', 'output/demo_HOSVD.avi');

% Perform mode-3 rank-1 partial svd
[core, U] = tensor_hosvd(T, 0, [0 0 1]);
L = tensor_ihosvd(core,U);
L = double(L);
S = double(T) - L;
