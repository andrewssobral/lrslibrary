function sol = shrinkage( p, alpha )
% For L1 regularization
% solves    min_x 0.5*|x-p|^2 + alpha*|x|
% p is a tenmat

pv = double( p );
temp = abs(pv) - alpha;
temp( temp < 0 ) = 0;
solv = temp .* sign(pv);
sol = tenmat( solv, p.rdims, p.cdims, p.tsize );

end