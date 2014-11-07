function Tv = shrinkage_v( Xv, alpha )
% Xv is a vector
% Tv is a vector

temp = abs(Xv) - alpha;
temp( temp < 0 ) = 0;
Tv = temp .* sign(Xv);

end