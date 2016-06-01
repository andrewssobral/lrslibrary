% TD | CP2 | PARAFAC2 decomposition
% process_video('TD', 'CP2', 'dataset/demo.avi', 'output/demo_CP2.avi');

r = 10;
A = double(T);
[B,H,C,P] = parafac2(A,r,[],[0 0 0 0 1]);
%%% PARAFAC2 reconstruction
for i = 1:size(C,1)
  L(:,:,i) = B*diag(C(i,:))*(P{i}*H)';
end
S = A - L;
