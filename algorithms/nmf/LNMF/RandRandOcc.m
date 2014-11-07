function [fea_o] = RandRandOcc(fea, h_img, w_img, ocld)

% random position occluding
% author@Guan Naiyang

n = size(fea,1);
for i = 1:n
    a = ceil(rand*(h_img+1-ocld));
    b = ceil(rand*(w_img+1-ocld));
    index = repmat(1:ocld, ocld, 1)'+repmat(0:h_img:(ocld-1)*h_img, ocld, 1);
    index = (b-1)*h_img+a-1+reshape(index, 1, ocld*ocld);
    fea(i, index) = 0;
end
fea_o = fea;

return;