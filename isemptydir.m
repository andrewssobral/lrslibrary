function x = isemptydir(p)
  if(isdir(p))
    f = dir(p);
    x = ~(length(f) > 3);
  else
    error('Error: %s is not a directory');
  end
end
