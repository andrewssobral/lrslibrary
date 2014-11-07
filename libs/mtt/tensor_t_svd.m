%% [U,S,V] = tensor_t_svd(A)
% 
function [U,S,V] = tensor_t_svd(A)
  D = fft(A,[],3);
  %D = A;
  n3 = size(A,3);
  
  for i = 1:n3 %nframes
    [Ux,Sx,Vx] = svd(D(:,:,i));
    Uy(:,:,i) = Ux;
    Sy(:,:,i) = Sx;
    Vy(:,:,i) = Vx;
  end
  
  U = ifft(Uy,[],3);
  S = ifft(Sy,[],3);
  V = ifft(Vy,[],3);
  %U = Uy;
  %S = Sy;
  %V = Vy;
end
