function spy3(T,varargin)
%SPY3 Visualize a third-order tensor's sparsity pattern.
%   spy3(T) plots the sparsity pattern of the third-order tensor T.
%
%   spy3(T,varargin) passes the parameters varargin{:} to the plot.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Check the dimensions.
if ndims(T) <= 2
    spy(T,varargin{:}); return;
elseif ndims(T) ~= 3
    error('spy3:ndims','ndims(T) should be 3.');
end

% Get the current axes.
cax = newplot;

% Set marker size.
units = get(gca,'units');
set(gca,'units','points');
pos = get(gca,'position');
ms = max(4,min(14,round(6*min(pos(3:4))/max(size(T)))));
set(gca,'units',units);

% Plot nonzero elements of T.
ind = find(T);
[i,j,k] = ind2sub(size(T),ind);
plot3(k,j,i,'marker','.','markersize',ms,'linestyle','none',varargin{:});

% Set axis properties.
xlabel('k');
ylabel('j');
zlabel('i');
xlim([1-5*eps size(T,3)]);
ylim([1-5*eps size(T,2)]);
zlim([1-5*eps size(T,1)]);
set(cax,'YDir','reverse');
set(cax,'ZDir','reverse');

% Display grid.
step = max(1,round(size(T)/6));
set(cax,'XTickMode','manual','YTickMode','manual','ZTickMode','manual');
set(cax,'XTick',[1:step(3):size(T,3)-1 size(T,3)]);
set(cax,'YTick',[1:step(2):size(T,2)-1 size(T,2)]);
set(cax,'ZTick',[1:step(1):size(T,1)-1 size(T,1)]);
grid on;

% Set title.
title(sprintf('nz = %i (%g%%)',nnz(T),nnz(T)/numel(T)*100));
