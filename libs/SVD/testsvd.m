%%
X = rand(100,400);
k = 100;

tic; [U,S,V] = svdsecon(X,k); toc; norm(U*S*V' - X,'fro')
tic; [U,S,V] = svds(X,k); toc; norm(U*S*V' - X,'fro')

tic; [U,S,V] = svdecon(X); toc; norm(U*S*V' - X,'fro')
tic; [U,S,V] = svd(X,'econ'); toc; norm(U*S*V' - X,'fro')

%%
%X = rand(400,100);
X = X';
k = 100;

tic; [U,S,V] = svdsecon(X,k); toc; norm(U*S*V' - X,'fro')
tic; [U,S,V] = svds(X,k); toc; norm(U*S*V' - X,'fro')

tic; [U,S,V] = svdecon(X); toc; norm(U*S*V' - X,'fro')
tic; [U,S,V] = svd(X,'econ'); toc; norm(U*S*V' - X,'fro')

%%
%X = rand(400,100);
X = X';
k = 9;

Xm = bsxfun(@minus,X,mean(X,2));
tic; [U,T] = pcaecon(X,k); toc
norm(Xm-U*T,'fro')

tic; [U,T] = pcasecon(X,k); toc
norm(Xm-U*T,'fro')

tic; [U,S,V] = svdsecon(Xm,k); toc
norm(Xm-U*S*V','fro')

tic; [U,S,V] = svdecon(Xm); toc
norm(Xm-U(:,1:k)*S(1:k,1:k)*V(:,1:k)','fro')

tic; [U,S,V] = svds(Xm,k); toc
norm(Xm-U*S*V','fro')

tic; [U,S,V] = svd(Xm); toc
norm(Xm-U(:,1:k)*S(1:k,1:k)*V(:,1:k)','fro')

