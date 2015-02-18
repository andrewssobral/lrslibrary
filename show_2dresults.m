%%
% 
function show_2dresults(m,n,varargin)
  k = nargin;
  vars = varargin;
  if(k < 3)
    error('Not supported operation!');
  end
  for i = 1:size(vars{1},2)
    for j = 1:k-2
      T = vars{j};
      T = reshape(T(:,i),m,n);
      switch(k)
        case 1
          imshow(T,[],'InitialMagnification','fit');
        case {2,3}
          subplot(1,k,j), imshow(T,[],'InitialMagnification','fit');
        case 4
          subplot(2,2,j), imshow(T,[],'InitialMagnification','fit');
        case {5,6}
          subplot(2,3,j), imshow(T,[],'InitialMagnification','fit');
        case {8}
          subplot(2,4,j), imshow(T,[],'InitialMagnification','fit');
        case {7,9}
          subplot(3,3,j), imshow(T,[],'InitialMagnification','fit');
        otherwise
          error('Not supported operation!');
      end
    end
    pause(0.01);
  end
end
