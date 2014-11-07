%%% void show_4dvideo(4dmatrix)
%
function show_4dvideo(video)
  for i = 1 : size(video,4)
    frame = video(:,:,:,i);
    imshow(frame,[],'InitialMagnification','fit');
    disp(i);
    pause(0.01);
  end
end
