% GA (Hauberg et al. 2014)
% process_video('RPCA', 'GA', 'dataset/demo.avi', 'output/demo_GA.avi');
L = grassmann_average(M', 1);
L = nma_rescale(L,min(M(:)),max(M(:))); % 0.968 0.207
L = repmat(L,1,size(M,2));
S = M - L;
% show_2dvideo(M,m,n);
% show_2dvideo(L,m,n);
% show_2dvideo(S,m,n);
