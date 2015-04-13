function slice3(T,varargin)
%SLICE3 Visualize a third-order tensor with slices.
%   slice3(T) visualizes the third-order tensor T by drawing its mode-1,
%   -2, and -3 slices using sliders to define their respective indices.
%   Press 'h' to show/hide the figure's controls.
%   
%   slice3(T,varargin) passes the parameters varargin{:} to the plot.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Check the dimensions.
if ndims(T) <= 2
    imagesc(T,varargin{:}); return;
elseif ndims(T) ~= 3
    error('slice3:T','ndims(T) should be 3.');
end

% If the data is complex, convert it to the modulus.
% Remove any NaNs and Infs.
if any(~isreal(T(:))), T = abs(T); end
T(isnan(T) | (isinf(T) & T < 0)) = min(T(~isinf(T(:))));
T(isinf(T) & T > 0) = max(T(~isinf(T(:))));

% Set up the figure.
idx = size(T); idx(2) = 1;
T = double(permute(T,[2 3 1]));
cax = newplot;
set(gcf,'Toolbar','figure');
set(datacursormode(gcf),'UpdateFcn',@datacursor);
zoom off; pan off; rotate3d off; datacursormode off;
obj{1} = uicontrol('Style','text','Position',[5 45 15 15],'String','i');
obj{2} = uicontrol('Style','slider','Position',[25 45 120 15], ...
                   'Min',1,'Max',size(T,3),'Value',size(T,3), ...
                   'SliderStep',1/(size(T,3)-1)*[1 1], ...
                   'Callback',{@redraw,1}, ...
                   'KeyPressFcn',@(obj,evt)toggle(evt.Key));
obj{3} = uicontrol('Style','text','Position',[5 25 15 15],'String','j');
obj{4} = uicontrol('Style','slider','Position',[25 25 120 15], ...
                   'Min',1,'Max',size(T,1),'Value',1, ...
                   'SliderStep',1/(size(T,1)-1)*[1 1], ...
                   'Callback',{@redraw,2}, ...
                   'KeyPressFcn',@(obj,evt)toggle(evt.Key));
obj{5} = uicontrol('Style','text','Position',[5 5 15 15],'String','k');
obj{6} = uicontrol('Style','slider','Position',[25 5 120 15], ...
                   'Min',1,'Max',size(T,2),'Value',size(T,2), ...
                   'SliderStep',1/(size(T,2)-1)*[1 1], ...
                   'Callback',{@redraw,3}, ...
                   'KeyPressFcn',@(obj,evt)toggle(evt.Key));
isVisible = true;
set(gcf,'KeyPressFcn',@(obj,evt)toggle(evt.Key));
redraw();

function toggle(key)
    if strcmpi(key,'h')
        if isVisible, val = 'off'; else val = 'on'; end
        for i = 1:length(obj), set(obj{i},'Visible',val); end
        isVisible = ~isVisible;
    end
end

function txt = datacursor(~,event_obj)
    pos = get(event_obj,'Position');
    txt = {['I: ' int2str(pos(3))], ...
           ['J: ' int2str(pos(2))], ...
           ['K: ' int2str(pos(1))], ...
           ['Value: ' num2str(T(pos(2),pos(1),pos(3)))]};
end

function redraw(hobj,~,ax)
    
    % Draw slices.
    if nargin >= 1, idx(ax) = round(get(hobj,'Value')); end
    slice(T,idx(3),idx(2),idx(1),varargin{:});
    shading flat;
    
    % Set axis properties.
    xlabel('k');
    ylabel('j');
    zlabel('i');
    xlim([1 size(T,2)]);
    ylim([1 size(T,1)]);
    zlim([1 size(T,3)]);
    set(cax,'YDir','reverse');
    set(cax,'ZDir','reverse');
    
    % Display grid.
    step = max(1,round(size(T)/6));
    set(cax,'XTickMode','manual','YTickMode','manual','ZTickMode','manual');
    set(cax,'XTick',[1:step(2):size(T,2)-1 size(T,2)]);
    set(cax,'YTick',[1:step(1):size(T,1)-1 size(T,1)]);
    set(cax,'ZTick',[1:step(3):size(T,3)-1 size(T,3)]);
    grid on;
    
end

end
