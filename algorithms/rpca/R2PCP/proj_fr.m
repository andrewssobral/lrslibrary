%% projection onto T_X(M_r)
%   X=U*diag(S)*V'
%   cf. [Van12]

function P=proj_fr(Z,U,V,sp)

if nargin<4, sp='t'; end

switch sp
    case 't'
        UtZ=U'*Z;
        P=U*UtZ+(Z*V)*V'-(U*(UtZ*V))*V';
        
    case 'n'
        UtZ=U'*Z;
        P=Z-(U*UtZ+(Z*V)*V'-(U*(UtZ*V))*V');
end


%% old version

% Pu=U*U';
% Pv=V*V';
% 
% switch sp
%     case 't'
%         P=Pu*Z+Z*Pv-Pu*Z*Pv;
%         
%     case 'n'
%         P=Z-Pu*Z-Z*Pv+Pu*Z*Pv;
% end
