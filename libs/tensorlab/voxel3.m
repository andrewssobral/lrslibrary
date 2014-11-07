function voxel3(T,varargin)
%VOXEL3 Visualize a third-order tensor with voxels.
%   voxel3(T) visualizes the third-order tensor T by plotting its elements
%   as voxels whose color and opacity are proportional to their value. The
%   figure contains two sliders for setting the parameters thresh and
%   degree (press 'h' to hide/show them). Let alpha = (T(i,j,k)-min(T(:))/
%   (max(T(:))-min(T(:))), then the opaqueness of each voxel, where 0 is
%   transparent and 1 is opaque, is computed as
%      
%      alpha                            if alpha >= thresh
%      thresh^(1-degree)*alpha^degree   if alpha <  thresh
%
%   voxel3(T,varargin) passes the parameters varargin{:} to the plot.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Check the dimensions.
if ndims(T) ~= 3
    error('voxel3:T','ndims(T) should be 3.');
end

% Check options.
idx = find(strcmpi('thresh',varargin),1,'last');
if idx < length(varargin)
    thresh = varargin{idx+1};
    varargin = [varargin(1:idx-1) varargin(idx+2:end)];
end
idx = find(strcmpi('degree',varargin),1,'last');
if idx < length(varargin)
    degree = varargin{idx+1};
    varargin = [varargin(1:idx-1) varargin(idx+2:end)];
end

% If the data is complex, convert it to the modulus.
% Remove any NaNs and Infs.
if any(~isreal(T(:))), T = abs(T); end
T(isnan(T) | (isinf(T) & T < 0)) = min(T(~isinf(T(:))));
T(isinf(T) & T > 0) = max(T(~isinf(T(:))));

% Compute the vertices.
st = size(T);
sv = st+1;
verts = zeros(prod(sv),3,'single');
verts(:,3) = repmat((.5:1:st(1)+.5).',size(verts,1)/sv(1),1);
verts(:,2) = repmat(kron((.5:1:st(2)+.5).',ones(sv(1),1)),sv(3),1);
verts(:,1) = kron((.5:1:st(3)+.5).',ones(size(verts,1)/sv(3),1));

% Compute the minimal number of colors and faces.
idx = {(1:prod(st+[1 0 0]))',(1:prod(st+[0 1 0]))',(1:prod(st+[0 0 1]))'};
off = [0 cumsum(cellfun(@numel,idx))];
color = zeros(sum(cellfun(@numel,idx)),1,'single');
faces = zeros(sum(cellfun(@numel,idx)),4,'single');
Tint = cat(1,T(1,:,:),max(T(1:end-1,:,:),T(2:end,:,:)),T(end,:,:));
[i,j,k] = ind2sub(size(Tint),idx{1});
color(off(1)+idx{1}) = single(Tint(:));
faces(off(1)+idx{1},:) = [sub2ind(sv,i,j  ,k)   sub2ind(sv,i,j+1,k) ...
                          sub2ind(sv,i,j+1,k+1) sub2ind(sv,i,j  ,k+1)];
Tint = cat(2,T(:,1,:),max(T(:,1:end-1,:),T(:,2:end,:)),T(:,end,:));
[i,j,k] = ind2sub(size(Tint),idx{2});
color(off(2)+idx{2}) = single(Tint(:));
faces(off(2)+idx{2},:) = [sub2ind(sv,i  ,j,k)   sub2ind(sv,i+1,j,k) ...
                          sub2ind(sv,i+1,j,k+1) sub2ind(sv,i  ,j,k+1)];
Tint = cat(3,T(:,:,1),max(T(:,:,1:end-1),T(:,:,2:end)),T(:,:,end));
[i,j,k] = ind2sub(size(Tint),idx{3});
color(off(3)+idx{3}) = single(Tint(:));
faces(off(3)+idx{3},:) = [sub2ind(sv,i  ,j  ,k) sub2ind(sv,i+1,j  ,k) ...
                          sub2ind(sv,i+1,j+1,k) sub2ind(sv,i  ,j+1,k)];

% Compute the opacity for each voxel.
mn = min(T(:));
mx = max(T(:));
alpha = (color-mn)/(mx-mn);
cutoff = 1e-4;

% Set up the figure.
cax = newplot;
if exist('degree','var') || exist('thresh','var')
    if exist('degree','var'), k = degree; else k = 1.0; end
    if exist('thresh','var'), t = thresh; else t = 0.5; end
else
    k = 1.0; t = 0.5;
    set(gcf,'Toolbar','figure');
    zoom off; pan off; rotate3d off; datacursormode off;
    lbl1 = uicontrol('Style','text','Position',[5 25 65 15], ...
                     'HorizontalAlignment','left');
    sld1 = uicontrol('Style','slider','Position',[75 25 120 15], ...
                     'Min',0,'Max',1,'Value',0.5, ...
                     'SliderStep',[0.1 0.25], ...
                     'Callback',{@redraw,'t'}, ...
                     'KeyPressFcn',@(obj,evt)toggle(evt.Key));
    lbl2 = uicontrol('Style','text','Position',[5 5 65 15], ...
                     'HorizontalAlignment','left');
    sld2 = uicontrol('Style','slider','Position',[75 5 120 15], ...
                     'Min',1,'Max',5,'Value',1, ...
                     'SliderStep',[0.25 0.5], ...
                     'Callback',{@redraw,'k'}, ...
                     'KeyPressFcn',@(obj,evt)toggle(evt.Key));
    isVisible = true;
    set(gcf,'KeyPressFcn',@(obj,evt)toggle(evt.Key));
end
redraw();

function toggle(key)
    % Toggle show or hide of uicontrols.
    if strcmpi(key,'h')
        onoff = 'on'; if isVisible, onoff = 'off'; end
        cellfun(@(c)set(c,'Visible',onoff),{lbl1,sld1,lbl2,sld2});
        isVisible = ~isVisible;
    end
end

function redraw(hobj,~,param)
    
    % Update the parameter.
    if nargin > 1
        switch param
            case 't', t = get(hobj,'Value');
            case 'k', k = get(hobj,'Value');
        end
    end
    if exist('lbl1','var')
        set(lbl1,'String',sprintf('thresh = %g',t));
        set(lbl2,'String',sprintf('degree = %g',k));
    end

    % Speed up rendering a bit.
    plot3(cax,nan,nan,nan,varargin{:});
    set(gcf,'DoubleBuffer','off');
    set(cax,'XLimMode','manual','YLimMode','manual', ...
            'ZLimMode','manual','CLimMode','manual','ALimMode','manual');
	
    % Set axis properties.
    xlabel('k');
    ylabel('j');
    zlabel('i');
    xlim([.5 st(3)+.5]);
    ylim([.5 st(2)+.5]);
    zlim([.5 st(1)+.5]);
    set(cax,'YDir','reverse');
    set(cax,'ZDir','reverse');
    caxis([mn mx]);

    % Display grid.
    step = max(1,round(st/6));
    set(cax,'XTickMode','manual','YTickMode','manual','ZTickMode','manual');
    tk = [1:step(3):st(3)-1 st(3)];
    if tk(end)-tk(end-1) < .3*step(3), tk = [tk(1:end-2) tk(end)]; end
    set(cax,'XTick',tk);
    tk = [1:step(2):st(2)-1 st(2)];
    if tk(end)-tk(end-1) < .3*step(2), tk = [tk(1:end-2) tk(end)]; end
    set(cax,'YTick',tk);
    tk = [1:step(1):st(1)-1 st(1)];
    if tk(end)-tk(end-1) < .3*step(1), tk = [tk(1:end-2) tk(end)]; end
    set(cax,'ZTick',tk);
    grid on;
    
    % Draw the voxels.
    opaq = alpha;
    opaq(alpha<t) = t^(1-k)*alpha(alpha<t).^k;
    opaq(opaq > 0.9) = 1; % Fix Matlab render bug.
    patch('Vertices',verts, ...
          'Faces',faces(opaq > cutoff,:), ...
          'EdgeColor','none','FaceAlpha','flat','FaceColor','flat', ...
          'FaceVertexAlphaData',opaq(opaq > cutoff), ...
          'FaceVertexCData',color(opaq > cutoff));
    
end

end
