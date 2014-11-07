load hall1-200;
matlabpool open;
X = normalize(XO);
[P Q] = RPMF(X, 2, 1, 1, 1e-2);
show(X, P * Q, abs(X - P * Q), [144 176]);
matlabpool close;
%%
L = P * Q;
%S = abs(X - P * Q);
S = X - P * Q;
isize = [144 176];
for i=1:min([300,size(X,2)])
    subplot(1,3,1);imshow(reshape(X(:,i),isize),[]), title('X(Sample)');
    subplot(1,3,2);imshow(reshape(L(:,i),isize),[]), title('L(Low-rank)');
    subplot(1,3,3);imshow(reshape(S(:,i),isize),[]), title('S(Sparse)');
    pause(0.01);
end