function data = gen_lowrank_tensor( Rs, rs )
% Rs and rs are vectors of ranks

% generate the core
N = length(Rs);
G = tensor( randn( rs ) );
U = cell( 1, N );

% generate orthonormal factors
for i = 1:N
    Ui = qmult( Rs(i) );
    U{i} = Ui(:,1:rs(i));  % Ri x ri
end

% multiply the factors to the core
data.X = ttm( G, U, 1:N );
data.rank = rs;

end