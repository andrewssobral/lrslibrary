%%% [2d_matrix, int, int, int] = convert_video3d_to_2d(3d_matrix)
%
function [M,m,n,p] = convert_video3d_to_2d(V)
  [m,n,p] = size(V);
  
  if(isa(V,'uint8'))
    M = uint8(zeros(m*n,p));
  else
    M = zeros(m*n,p);
  end
  
  for i = 1:p
    M(:,i) = reshape(V(:,:,i),[],1);
  end
end
