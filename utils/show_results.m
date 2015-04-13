%%% void show_results(matrix,matrix,matrix,matrix,int,int,int)
% I - Input sequence
% L - Low-rank sequence
% S - Sparse sequence
% O - Outliers sequence
% nFrame - Number of frames
% vidHeight - Video height
% vidWidth - Video width
%
function show_results(I,L,S,O,nFrame,vidHeight,vidWidth)
  warning('off','all');
  clf;
  for i = 1:nFrame
    Input = reshape(I(:,i),vidHeight,vidWidth);
    %Input = uint8(Input);
    
    if(size(L) == size(I))
      LowRank = reshape(L(:,i),vidHeight,vidWidth);
      %LowRank1 = reshape(L(:,i),vidHeight,vidWidth);
      %LowRank = mat2gray(LowRank1);
      %LowRank = im2uint8(LowRank);
    else
      LowRank = uint8(zeros(size(Input)));
    end

    if(size(S) == size(I))
      Sparse = reshape(S(:,i),vidHeight,vidWidth);
      %Sparse1 = reshape(S(:,i),vidHeight,vidWidth);
      %Sparse = mat2gray(Sparse1);
      %Sparse = im2uint8(Sparse);
    else
      Sparse = uint8(zeros(size(Input)));
    end

    if(~isempty(O))
      Outlier = reshape(O(:,i),vidHeight,vidWidth);
      %Outlier1 = reshape(O(:,i),vidHeight,vidWidth);
      %Outlier = mat2gray(Outlier1);
      %Outlier = im2uint8(Outlier);
      %Outlier = medfilt2(Outlier, [5 5]);
    else
      Outlier = uint8(zeros(size(Input)));
    end

    subplot(2,3,1), imshow(Input,[]), title('Input (I)');
    subplot(2,3,2), imshow(LowRank,[]), title('Low Rank (L)');
    subplot(2,3,3), imshow(Sparse,[]), title('Sparse (S)');
    subplot(2,3,4), imshow(Outlier,[]), title('Outliers (O)');
    subplot(2,3,5), imshow(medfilt2(Outlier, [5 5]),[]), title('Filtered Outliers');
    %disp(i);
    pause(0.01);
  end
  warning('on','all');
end

