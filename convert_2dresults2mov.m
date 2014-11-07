%%% movobj = convert_2dresults2mov(2dmatrix,2dmatrix,2dmatrix,2dmatrix,struct)
% I - Input sequence
% L - Low-rank sequence
% S - Sparse sequence
% O - Outliers sequence
%
function [movobj] = convert_2dresults2mov(I,L,S,O,video)
  %warning('off','all');
  displog('Converting results to movie...');
  frames = video.nrFramesTotal;
  height = video.height;
  width = video.width;
  
  if(~isempty(I)) movobj.I(1:frames) = struct('cdata', [], 'colormap', []); else movobj.I = []; end
  if(~isempty(L)) movobj.L(1:frames) = struct('cdata', [], 'colormap', []); else movobj.L = []; end
  if(~isempty(S)) movobj.S(1:frames) = struct('cdata', [], 'colormap', []); else movobj.S = []; end
  if(~isempty(O)) movobj.O(1:frames) = struct('cdata', [], 'colormap', []); else movobj.O = []; end
  
  for i = 1 : frames
    %
    % Convert input
    %
    if(~isempty(I))
      Input = reshape(I(:,i),height,width);
      Input = mat2gray(Input);
      Input = im2uint8(Input);
      rgbImage = repmat(Input,[1 1 3]);
      movobj.I(i).cdata = rgbImage;
    end
    %
    % Convert low-rank results
    %
    if(~isempty(L))
      LowRank = reshape(L(:,i),height,width);
      LowRank = mat2gray(LowRank);
      LowRank = im2uint8(LowRank);
      rgbImage = repmat(LowRank,[1 1 3]);
      movobj.L(i).cdata = rgbImage;
    end
    %
    % Store sparse results
    %
    if(~isempty(S))
      Sparse = reshape(S(:,i),height,width);
      Sparse = mat2gray(Sparse);
      Sparse = im2uint8(Sparse);
      rgbImage = repmat(Sparse,[1 1 3]);
      movobj.S(i).cdata = rgbImage;
    end
    %
    % Store foreground results
    %
    if(~isempty(O))
      Outlier = reshape(O(:,i),height,width);
      Outlier = mat2gray(Outlier);
      Outlier = im2uint8(Outlier);
      Outlier = medfilt2(Outlier, [5 5]);
      rgbImage = repmat(Outlier,[1 1 3]);
      movobj.O(i).cdata = rgbImage;
    end
  end
  
  %warning('on','all');
end
