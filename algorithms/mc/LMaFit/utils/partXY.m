function v = partXY(Xt,Y,I,J,L)

% Input:
%          Xt --- k x m matrix
%           Y --- k x n matrix
%           I --- row indices
%           J --- column indices
%           L --- length of v
% Output:
%           v --- XY(I,J);

v = zeros(L,1);

for p = 1:L
    i = I(p);
    j = J(p);
    v(p) = Xt(:,i)'*Y(:,j);    
end
