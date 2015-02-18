addpath ./propack
addpath ./yima_rpca
addpath ./lib

X = rand(1000, 1000);
k = 20;

tic;[U1 s1 V1] = svdex(X, k, struct('svd_solver', 'svd'));t1 = toc;
tic;[U2 s2 V2] = svdex(X, k, struct('svd_solver', 'svds'));t2 = toc;
tic;[U3 s3 V3] = svdex(X, k, struct('svd_solver', 'lansvd', 'init', randn(size(X,1),1)));t3 = toc;
figure; plot(1:k,s1,'r', 1:k,s2,'g^', 1:k,s3,'bv');
fprintf('Time:\n     svd=%fs\n    svds=%fs\n  lansvd=%fs\n',t1,t2,t3);

Y1 = bsxfun(@times, U1, s1')*V1';
err1 = max(abs(X(:) - Y1(:)));
Y2 = bsxfun(@times, U2, s2')*V2';
err2 = max(abs(X(:) - Y2(:)));
Y3 = bsxfun(@times, U3, s3')*V3';
err3 = max(abs(X(:) - Y3(:)));
fprintf('Error:\n     svd=%f\n    svds=%f\n  lansvd=%f\n',err1,err2,err3);
