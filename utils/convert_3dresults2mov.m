%%% movobj = convert_3dresults2mov(3dmatrix,3dmatrix,3dmatrix,3dmatrix,struct)
% I - Input sequence
% L - Low-rank sequence
% S - Sparse sequence
% O - Outliers sequence
%
function [movobj] = convert_3dresults2mov(I,L,S,O,nframes)
  %warning('off','all');
  displog('Converting results to movie...');
  
  if(~isempty(I)) movobj.I(1:nframes) = struct('cdata', [], 'colormap', []); else movobj.I = []; end
  if(~isempty(L)) movobj.L(1:nframes) = struct('cdata', [], 'colormap', []); else movobj.L = []; end
  if(~isempty(S)) movobj.S(1:nframes) = struct('cdata', [], 'colormap', []); else movobj.S = []; end
  if(~isempty(O)) movobj.O(1:nframes) = struct('cdata', [], 'colormap', []); else movobj.O = []; end
  
  for i = 1 : nframes
    %
    % Convert input
    %
    if(~isempty(I))
      Input = I(:,:,i);
      Input = mat2gray(Input);
      Input = im2uint8(Input);
      rgbImage = repmat(Input,[1 1 3]);
      movobj.I(i).cdata = rgbImage;
    end
    %
    % Convert low-rank results
    %
    if(~isempty(L))
      LowRank = L(:,:,i);
      LowRank = mat2gray(LowRank);
      LowRank = im2uint8(LowRank);
      rgbImage = repmat(LowRank,[1 1 3]);
      movobj.L(i).cdata = rgbImage;
    end
    %
    % Store sparse results
    %
    if(~isempty(S))
      Sparse = S(:,:,i);
      Sparse = mat2gray(Sparse);
      Sparse = im2uint8(Sparse);
      rgbImage = repmat(Sparse,[1 1 3]);
      movobj.S(i).cdata = rgbImage;
    end
    %
    % Store foreground results
    %
    if(~isempty(O))
      Outlier = O(:,:,i);
      Outlier = mat2gray(Outlier);
      Outlier = im2uint8(Outlier);
      Outlier = medfilt2(Outlier, [5 5]);
      rgbImage = repmat(Outlier,[1 1 3]);
      movobj.O(i).cdata = rgbImage;
    end
  end
  
  %warning('on','all');
end
