function show(X, L, S, isize)
    figure;
    for i=1:min([300,size(X,2)])
        subplot(1,3,1);imshow(reshape(X(:,i),isize),[]);%title('X(Sample)');
        subplot(1,3,2);imshow(reshape(L(:,i),isize),[]);%title('L(Low-rank)');
        subplot(1,3,3);imshow(reshape(S(:,i),isize),[]);%title('S(Sparse)');
        pause(0.01);
    end
end