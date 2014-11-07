%% [core,U] = tensor_hosvd(T,t,r)
%
%%% Normal HoSVD
% [core,U] = tensor_hosvd(T)
% [core,U] = tensor_hosvd(T,0)
% [core,U] = tensor_hosvd(T,[0 0 0])
%
%%% Truncate mode-3 basis matrice (keep the first 2 eigenvalues)
% [core,U] = tensor_hosvd(T,[0 0 2])
%
%%% Perform mode-3 rank-2 partial svd
% [core,U] = tensor_hosvd(T,0,[0 0 2])
%
function [core,U] = tensor_hosvd(T,t,r)
  core = T;
  
  if(nargin < 3) r = zeros(1,ndims(T)); end
  if(nargin < 2) t = zeros(1,ndims(T)); end
  
  if(t == 0) t = zeros(1,ndims(T)); end
  if(r == 0) r = zeros(1,ndims(T)); end
  
  for i = 1:ndims(T)
    M{i} = double(tenmat(T,i));
    %M{i} = unfolding(double(T),i);
    
    %%% perform partial svd
    if(r(i) > 0)
      [U{i},S{i},V{i}] = svds(M{i},r(i),'L');
    else
      [U{i},S{i},V{i}] = svd(M{i});
    end

    %%% truncate basis matrice
    if(t(i) > 0)
      U{i}(:,t(i)+1:end) = 0;
    end

    %[m n] = size(U{i});
    %s(i) = n;
    %core = tensor(folding(U{i}'*unfolding(double(core),i), i, s));
    %core = double(ttm(tensor(core),U{i},i,'t'));
    core = ttm(core,U{i}',i);
  end
end
