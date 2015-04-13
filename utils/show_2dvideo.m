%%% show_2dvideo(2dmatrix,int,int)
%
function show_2dvideo(M,m,n)
  for i = 1 : size(M,2)
    I = reshape(M(:,i),m,n);
    imshow(I,[],'InitialMagnification','fit');
    disp(i);
    pause(0.01);
  end
end
