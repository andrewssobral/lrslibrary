%%% void convert_video4d_to_avi(V,filename)
%
function convert_video4d_to_avi(video,filename)
  movobj(1:size(video,4)) = struct('cdata',[], 'colormap',[]);
  
  for i = 1:size(video,4)
    frame = video(:,:,:,i);
    
    if(size(frame,3) == 1)
      frame = repmat(frame,[1 1 3]);
    end
    
    movobj(i).cdata = frame;
  end
  
  movie2avi(movobj, filename, 'compression', 'None');
end
