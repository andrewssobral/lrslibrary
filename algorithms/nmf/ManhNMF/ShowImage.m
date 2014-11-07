function ShowImage(data, height, width, nRow, nColumn, reverse)

[p, m] = size(data);
if width * height ~= p
    error('<EigFace>: incorrect width or height.\n');
end

if nColumn * nRow > m
    error('<EigFace>: incorrect column number or row number.\n');
end

lineW = 1;
horiW = lineW;
vertW = lineW;
% scaling data
data = data/max(data(:));
canvH = height*nRow+horiW*(nRow+1);
canvW = width*nColumn+vertW*(nColumn+1);

% construct canvas
canvColor = .5;
if reverse
    canvColor = 1-canvColor;
end
Y = canvColor*ones(canvH,canvW);

% draw the datas
vertProb = lineW;
for i=1:nRow
    horiProb = lineW;
    for j=1:nColumn
        if reverse
            Y(vertProb+(1:height),  horiProb+(1:width)) = 1-reshape(data(:,(i-1)*nColumn+j),[height,width]);
        else
            Y(vertProb+(1:height),  horiProb+(1:width)) = reshape(data(:,(i-1)*nColumn+j),[height,width]);
        end
        % draw dashed line
        if j~=nColumn
            Y(vertProb+(1:floor(height/3):height), horiProb+width+1) = ~canvColor;
        end
        if i~=nRow
            Y(vertProb+height+1, horiProb+(1:floor(width/3):width)) = ~canvColor;
        end
        % move horizontal probe
        horiProb = horiProb+width+vertW;
    end
    % move vertical probe
    vertProb = vertProb+height+horiW;
end 

% rescaling
Y = Y*255;

% show the image
imagesc(Y);
colormap(gray);
axis image off;