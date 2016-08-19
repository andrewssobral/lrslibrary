function [bnd,gap] = refinebounds(D,bnd,tol1)
%REFINEBONDS  Refines error bounds for Ritz values based on gap-structure
% 
%  bnd = refinebounds(lambda,bnd,tol1) 
%
%  Treat eigenvalues closer than tol1 as a cluster.

% Rasmus Munk Larsen, DAIMI, 1998

j = length(D);

if j<=1
  return
end
% Sort eigenvalues to use interlacing theorem correctly
[D,PERM] = sort(D);
bnd = bnd(PERM);


% Massage error bounds for very close Ritz values
eps34 = sqrt(eps*sqrt(eps));
[y,mid] = max(bnd);
for l=[-1,1]    
  for i=((j+1)-l*(j-1))/2:l:mid-l
    if abs(D(i+l)-D(i)) < eps34*abs(D(i))
      if bnd(i)>tol1 & bnd(i+l)>tol1
	bnd(i+l) = pythag(bnd(i),bnd(i+l));
	bnd(i) = 0;
      end
    end
  end
end
% Refine error bounds
gap = inf*ones(1,j);
gap(1:j-1) = min([gap(1:j-1);[D(2:j)-bnd(2:j)-D(1:j-1)]']);
gap(2:j) = min([gap(2:j);[D(2:j)-D(1:j-1)-bnd(1:j-1)]']);
gap = gap(:);
I = find(gap>bnd);
bnd(I) = bnd(I).*(bnd(I)./gap(I));

bnd(PERM) =  bnd;