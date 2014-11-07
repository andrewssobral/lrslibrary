function [Iwarp,Omega] = warpImg(I,tau,mode,extrapval)

% This function warps I to get a new image Iwarp based on different
% transformation
% tau --- the transformation parameter
% length(tau) = 3: rigid
% length(tau) = 4: similarity
% length(tau) = 6: affine
% length(tau) = 8: projective
% mode = 0: backward [x';y'] = f([x;y];tau); I2(x,y) = I1(x',y');
% mode = 1: foreward [x';y'] = f([x;y];tau); I2(x',y') = I1(x,y);
% extrapval: a scaler to specify the extrapolation value
% if extrapval = [], the nearest value from the border is copied
% The origin of the field is at the center of the image!
% RIGID
%   c(tau(1))   -s(tau(1))  tau(2)
%   s(tau(1))   c(tau(1))   tau(3)
% SIMILARITY
%   1+tau(1)    -tau(2)     tau(3)
%   tau(2)      1+tau(1)    tau(4)
% AFFINE
%   1+tau(1)    tau(3)      tau(5)
%   tau(2)      1+tau(4)    tau(6)
% PROJECTIVE
%   1+tau(1)    tau(4)      tau(7)
%   tau(2)      1+tau(5)    tau(8)
%   tau(3)      tau(6)      1

if nargin < 3 
    mode = 0;
end
if nargin < 4 
    extrapval = [];
end
I = double(I);
sizeI = size(I);

% grids before transfomation
[Yb4,Xb4] = meshgrid(1:sizeI(2),1:sizeI(1));
% coordinates before transformation
Coordb4 = [ Xb4(:)'- round(sizeI(1)/2); ...
            Yb4(:)'- round(sizeI(2)/2); ...
            ones(1,length(Xb4(:)))]; % 3-by-L [x;y;1]
        
% coordinates after transformation
switch length(tau)
    case 3
        M = [cos(tau(1)) -sin(tau(1)) tau(2);...
             sin(tau(1)) cos(tau(1))  tau(3);...
             0 0 1];
        if mode == 0 
            Coordaft = M*Coordb4;
        else
            Coordaft = M\Coordb4;
        end
    case 4 
        M = [1+tau(1)    -tau(2)     tau(3);...
             tau(2)      1+tau(1)    tau(4)];
        if mode == 0 
            Coordaft = M*Coordb4;
        else
            Coordaft = M\Coordb4;
        end
    case 6
        M = [1+tau(1)    tau(3)      tau(5);...
             tau(2)      1+tau(4)    tau(6);...
             0 0 1];
        if mode == 0 
            Coordaft = M*Coordb4;
        else
            Coordaft = M\Coordb4;
        end
    case 8
        M = [1+tau(1)    tau(4)      tau(7);...
             tau(2)      1+tau(5)    tau(8);...
             tau(3)      tau(6)      1];
        if mode == 0 
            Coordaft = M*Coordb4;
        else
            Coordaft = M\Coordb4;
        end
        Coordaft = bsxfun(@rdivide,Coordaft,Coordaft(3,:)+eps);
    otherwise
        error('length(tau) is illegal');
end
% grids after transfomation
Xaft = reshape(Coordaft(1,:)',sizeI) + round(sizeI(1)/2);
Yaft = reshape(Coordaft(2,:)',sizeI) + round(sizeI(2)/2);

% if move out of the image plane
Omega = false(sizeI);
Omega(Xaft<1|Xaft>sizeI(1)|Yaft<1|Yaft>sizeI(2)) = true;

% interpolation
if isempty(extrapval)
    Xaft(Xaft<1) = 1;
    Xaft(Xaft>sizeI(1)) = sizeI(1);
    Yaft(Yaft<1) = 1;
    Yaft(Yaft>sizeI(2)) = sizeI(2);
    Iwarp = interp2(Yb4,Xb4,I,Yaft,Xaft,'*bicubic');
else
    Iwarp = interp2(Yb4,Xb4,I,Yaft,Xaft,'*bicubic',extrapval);
end