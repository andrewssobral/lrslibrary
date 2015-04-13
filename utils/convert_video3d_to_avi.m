%%% void convert_video3d_to_avi(V,filename)
%
function convert_video3d_to_avi(V,filename)
  movobj(1:size(V,3)) = struct('cdata',[], 'colormap',[]);
  for i = 1:size(V,3)
    frame = V(:,:,i);
    frame = mat2gray(frame);
    frame = im2uint8(frame);
    rgb = repmat(frame,[1 1 3]);
    movobj(i).cdata = rgb;
  end
  movie2avi(movobj, filename, 'compression', 'None');
end
