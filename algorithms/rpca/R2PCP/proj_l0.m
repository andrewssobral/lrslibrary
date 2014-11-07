function Y=proj_l0(Z,X,sp)

if nargin<3, sp='t'; end

Y=Z;
switch sp
    case 't'
        Y(X==0)=0;
    case 'n'
        Y(X~=0)=0;
end
