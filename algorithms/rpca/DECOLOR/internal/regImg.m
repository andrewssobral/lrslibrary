function [I2warp,tau,residue,OmegaOut] = regImg(I1,I2,tau,weight,maxIts)

% This function registers I2 to I1 based on different different transformation
% tau --- the transformation parameter,
% which must be initialized, since its length specifies the type of transf.
% length(tau) = 3: rigid
% length(tau) = 4: similarity
% length(tau) = 6: affine
% length(tau) = 8: projective

% initialize variables
I1 = double(I1);
I2 = double(I2);
sizeI = size(I1);
sizeD = sizeI(1)*sizeI(2);

if nargin < 3 || isempty(tau)
    tau = zeros(6,1); % the default is affine transf.
end
if  nargin < 4 || isempty(weight)
    weight = ones(sizeD,1);
else
    weight = reshape(weight,[sizeD,1]);
end
if  nargin < 5 || isempty(maxIts)
    maxIts = 50;
end

for iter = 1:maxIts
    
    % warping to calculate I2(X) = I1(f(X,tau))
    [I2warp,OmegaOut] = warpImg(I2,tau);
    % get derivatives
    [I2warp_x,I2warp_y] = getGradient(I2warp);
    
    %% compute least square solution
    y = reshape(I1-I2warp,[sizeD,1]);
    
    % coordinates
    [yCoord,xCoord] = meshgrid(1:sizeI(2),1:sizeI(1));
    xCoord = xCoord - round(sizeI(1)/2);
    yCoord = yCoord - round(sizeI(2)/2);
    
    % compute X(i,:) = [I2x,I2y]*Jacob(d[x;y]/d tau)
    switch length(tau)
        case 3
            X = [(-sin(tau(1))*xCoord(:)-cos(tau(1))*yCoord(:)).*I2warp_x(:)...
                +( cos(tau(1))*xCoord(:)-sin(tau(1))*yCoord(:)).*I2warp_y(:),...
                I2warp_x(:), I2warp_y(:)]; % 1-by-3
        case 4
            X = [xCoord(:).*I2warp_x(:)+yCoord(:).*I2warp_y(:),...
                -yCoord(:).*I2warp_x(:)+xCoord(:).*I2warp_y(:),...
                I2warp_x(:), I2warp_y(:)]; % 1-by-4
        case 6
            X = [xCoord(:).*I2warp_x(:), xCoord(:).*I2warp_y(:),...
                 yCoord(:).*I2warp_x(:), yCoord(:).*I2warp_y(:),...
                 I2warp_x(:), I2warp_y(:)]; % 1-by-6
        case 8
            X = [ xCoord(:).*I2warp_x(:), xCoord(:).*I2warp_y(:), -(xCoord(:).^2).*I2warp_x(:)-(xCoord(:).*yCoord(:)).*I2warp_y(:),...
                  yCoord(:).*I2warp_x(:), yCoord(:).*I2warp_y(:), -(yCoord(:).^2).*I2warp_y(:)-(xCoord(:).*yCoord(:)).*I2warp_x(:),...
                  I2warp_x(:), I2warp_y(:)]; % 1-by-6
            D = [xCoord(:),yCoord(:)]*[tau(3);tau(6)]+1;
            X = bsxfun(@rdivide,X,D+eps);
        otherwise
    end
    
    % solve
    A = X'*bsxfun(@times,X,weight);
    dtau = (A+0.001*diag(diag(A)))\(X'*(weight.*y));
    tau = tau + dtau;

    % check termination condition
    if max(abs(dtau./(tau+eps))) < 0.001;
        break
    end
end

residue = reshape( y - X*dtau, sizeI );
end





