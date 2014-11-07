%% TENSOR MATRICIZATION
% V1 - mode-1 matricization
% V2 - mode-2 matricization
% V3 - mode-3 matricization
%
function [V1,V2,V3] = tensor_matricization(V)
  %%% mode-1 matricization (#lins)
  % V1 -> R{dim1,dim2*dim3} = R{60,80x52} = R{60,4160}
  V1 = reshape(V,size(V,1),[]);
  % imshow(V1,[])

  %%% mode-2 matricization (#cols)
  % V2 -> R{dim2,dim1*dim3} = R{80,60x52} = R{60,3120}
  V2 = permute(V,[2 1 3]);
  V2 = reshape(V2,size(V,2),[]);
  % imshow(V2,[])

  %%% mode-3 matricization (#deph)
  % V3 -> R{dim3,dim1*dim2} = R{52,60x80} = R{60,4800}
  V3 = permute(V,[3 1 2]);
  V3 = reshape(V3,size(V,3),[]);
  % imshow(V3,[])
end

