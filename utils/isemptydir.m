function x = isemptydir(p,k)
  if(isdir(p))
    f = dir(p);
    x = ~(length(f) > k);
  else
    error('Error: %s is not a directory');
  end
end
