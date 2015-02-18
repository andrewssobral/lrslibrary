function [ ims ] = ShowImages( x, sz, cols )
%[ ims ] = showimages( x,sz,cols )
%SHOWIMAGES show images
%   x: images, each column as an image
%	sz: [height weight]
%	cols: nubmer of images in a row
% author: Liang Xiong (lxiong@cs.cmu.edu)

cla;
if isempty(x)
    return;
end

height = sz(1);
width = sz(2);

x = full(x);

if nargin < 3
    cols = min(size(x,2),floor(800/width));
end		
rows = ceil(size(x,2)/cols);

ims = zeros(height*rows,width*cols);
for ind = 1:size(x,2)
	yy = ceil(ind/cols);
	xx = ind - cols*(yy - 1);
	ims(((yy-1)*height+1):(yy*height),((xx-1)*width+1):(xx*width)) = ...
            reshape(x(:,ind),height,width);
end

mi = min(x(:));	
ma = max(x(:)) + 1e-2;
if nargout == 0
    imshow(ims,[mi ma]);
end
