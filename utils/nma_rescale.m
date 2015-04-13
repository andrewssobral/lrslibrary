function C = nma_rescale(A,new_min,new_max)
  %Nasser M. Abbasi 011212
  %NO ERROR CHECKING DONE ON INPUT. Rescale a matrix or a vector A
  current_max = max(A(:));
  current_min = min(A(:));
  C =((A-current_min)*(new_max-new_min))/(current_max-current_min) + new_min;
end