%%% void show_3dvideo(3dmatrix)
%
function show_3dvideo(V)
  for i = 1 : size(V,3)
    frame = V(:,:,i);
    imshow(frame,[],'InitialMagnification','fit');
    disp(i);
    pause(0.01);
  end
end
