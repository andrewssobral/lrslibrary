%% show_3dtensors(T1,...,Tn)
% 
function show_3dtensors(varargin)
  n = nargin;
  % display(['Number of tensors: ' num2str(n)]);
  V = varargin{1};
  for i = 1:size(V,3)
    for j = 1:n
      T = varargin{j};
      if(isa(T,'tensor'))
        T = double(T);
      end
      switch(n)
        case 1
          imshow(T(:,:,i),[],'InitialMagnification','fit');
        case {2,3}
          subplot(1,n,j), imshow(T(:,:,i),[],'InitialMagnification','fit');
        case 4
          subplot(2,2,j), imshow(T(:,:,i),[],'InitialMagnification','fit');
        case {5,6}
          subplot(2,3,j), imshow(T(:,:,i),[],'InitialMagnification','fit');
        case {8}
          subplot(2,4,j), imshow(T(:,:,i),[],'InitialMagnification','fit');
        case {7,9}
          subplot(3,3,j), imshow(T(:,:,i),[],'InitialMagnification','fit');
        otherwise
          error('Not supported operation!');
      end
    end
    pause(0.01);
  end
end
