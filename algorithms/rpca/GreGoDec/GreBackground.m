function [L,S,G,error,time]=GreBackground(X,rank,tau,power,isize,tol,k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Background Modeling by Greedy Semi-Soft GoDec (GreBsmo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%isize: image size of each frame, for example, 1024x768
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tianyi Zhou, 2013, All rights reserved.

[L,S,error,time]=GreGoDec(X,rank,tau,tol,power,k);
G=X-L-S;

figure;
for i=1:min([200,size(X,2)])
    subplot(2,2,1);imagesc(reshape(X(:,i),isize));colormap(gray);axis image;axis off;title('X(Sample)');
    subplot(2,2,2);imagesc(reshape(L(:,i),isize));colormap(gray);axis image;axis off;title('L(Low-rank)');
    %subplot(1,4,3);imagesc(reshape(X(:,i).*(~~S(:,i)),isize));colormap(gray);axis image;axis off;title('S(Sparse)');
    subplot(2,2,4);imagesc(reshape(S(:,i),isize));colormap(gray);axis image;axis off;title('S(Sparse)');
    subplot(2,2,3);imagesc(reshape(G(:,i),isize));colormap(gray);axis image;axis off;title('G(Noise)');
    pause(0.1);
    if i<10
        print('-djpeg', '-r250', sprintf('%s%s','00',num2str(i)));
    elseif i<100
        print('-djpeg', '-r250', sprintf('%s%s','0',num2str(i)));
    else
        print('-djpeg', '-r250', sprintf('%s',num2str(i)));
    end
end