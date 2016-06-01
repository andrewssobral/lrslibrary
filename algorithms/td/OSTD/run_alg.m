nframes = size(T,3);
for i = 1:nframes
  disp(['processing frame: ' num2str(i)]);
  if(ndims(T) == 3)
    frame = T(:,:,i);
  end
  if(ndims(T) == 4)
    frame = T(:,:,:,i);
  end
  T_i = tensor(frame);
  if(i == 1) Tm = []; end
  [Tlowrank,Tsparse,Tmask,Tm] = OSTD(T_i,i,Tm);
  L(:,:,i) = Tlowrank;
  S(:,:,i) = Tsparse;
  O(:,:,i) = Tmask;
end
