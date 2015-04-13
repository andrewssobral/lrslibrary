%%% void show_Svideo(struct)
%
function show_Svideo(V)
  for i = 1 : size(V,2)
    frame = V{i};
    imshow(frame,[],'InitialMagnification','fit');
    disp(i);
    pause(0.01);
  end
end
