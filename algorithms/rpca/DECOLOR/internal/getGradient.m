function [phi_x phi_y] = getGradient(phi,sigma)

% get the partial derivative of a 2D function phi
% Input:
% 	phi: 2D function
% Output:
%	phi_x:	partial derivative of x direction
%	phi_y:  partial derivative of y direction
%	if index is not specified, the output will be reshaped as the same dimension of the input phi
% Author: Zhou Xiaowei, ECE@HKUST, eexwzhou@ust.hk, 2010-12-21

if exist('sigma','var')
    phi = gSmooth(phi,sigma);
end

detector = [1 -8 0 8 -1]'/12;
phi_x = imfilter(phi,detector, 'replicate');
phi_y = imfilter(phi,detector','replicate');

end

function I = gSmooth(I,sigma)

dim = ndims(I);
if dim ~= length(sigma)
    error('the length of sigma should be the same with the dimension of I');
end
mask = getGaussian(size(I),sigma);
I = real(ifftshift(ifftn(fftn(mask).*fftn(I))));
end

function mask = getGaussian(size,sigma)
% This function generates a 2D Gaussian mask for filtering
% Input:
%   size:   a vecter giving the size of the mask along x and y direciton
%   sigma:  a vector defining the sigma value along x and y direction
% Output:
%   mask:   a 2D matrix whose enrty is the coefficient

if length(size) ~= length(sigma)
	error('the vector size and sigma should have the same length');
else
	dim = length(size);
end

sigma(sigma==0) = 0.1;

switch dim
	case 1
		X = [1:size]-(1+size)/2;
		mask = exp(-(X/(sigma+eps)).^2);
		mask = mask/sum(mask);
	case 2
		X = ([1:size(1)]'-(1+size(1))/2)*ones(1,size(2));
		Y = ones(size(1),1)*([1:size(2)]-(1+size(2))/2);
		mask = exp(-((X/(sigma(1)+eps)).^2+(Y/(sigma(2)+eps)).^2)/2 );
		mask = mask/sum(mask(:));
	case 3
		xind = [1:size(1)]'-(1+size(1))/2;
		X = zeros(size);
		for i = 1:length(xind)
			X(i,:,:) = xind(i);
		end
		yind = [1:size(2)]'-(1+size(2))/2;
		Y = zeros(size);
		for i = 1:length(yind)
			Y(:,i,:) = yind(i);
		end
		zind = [1:size(3)]'-(1+size(3))/2;
		Z = zeros(size);
		for i = 1:length(zind)
			Z(:,:,i) = zind(i);
		end
		mask = exp(-((X/(sigma(1)+eps)).^2+(Y/(sigma(2)+eps)).^2+(Z/(sigma(3)+eps)).^2)/2);
		mask = mask/sum(mask(:));
	otherwise
		error('can only support dimension 3 below');
end

end

    