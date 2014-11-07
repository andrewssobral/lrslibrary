function [size_core,L_ub,L_tmlsvd] = mlrankest(T,options)
%MLRANKEST Estimate multilinear rank.
%   mlrankest(T) plots an L-curve of Pareto optimal multilinear ranks for
%   low multilinear rank approximation, based on an upper bound of the
%   truncation error. The x-axis corresponds to the compression ratio
%   obtained by a LMLRA for a given multilinear rank, and the y-axis
%   corresponds to the relative error of that approximation. Additionally,
%   the corner of the resulting L-curve is estimated and the relative error
%   of a truncated MLSVD in the neighbourhood of the corner is computed.
%
%   size_core = mlrankest(T) does not plot anything and instead returns the
%   vector of truncation ranks corresponding to the corner of the L-curve.
%
%   [size_core,L_ub,L_tmlsvd] = mlrankest(T) also returns the set of Pareto
%   optimal multilinear ranks L_ub(:,1:ndims(T)) and the upper bound on
%   their truncation error L_ub(:,end), as well as the multilinear ranks
%   L_tmlsvd(:,1:ndims(T)) in a neighbourhood around the L-curve corner and
%   their corresponding relative errors of the truncated MLSVD
%   L_tmlsvd(:,end).
%
%   mlrankest(T,options) may be used to set the following options:
%
%      options.Range = 10      - The range of truncation ranks near the
%                                L-curve corner of which to compute the
%                                truncated MLSVD error. Set a negative
%                                range to skip this step.
%      options.XMultiplier = 1 - The importance of the compression ratio in
%                                determining the L-curve corner, relative
%                                to the importance of the relative error of
%                                the approximation.
%
%   See also rankest.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] J.L. Castellanos, S. Gomez, V. Guerra, "The triangle method for
%       finding the corner of the L-curve," Applied Numerical Mathematics,
%       Vol. 43, No. 4, 2002, pp. 359-373.

% Incomplete/sparse tensors not supported yet.
T = fmt(T,true);
if isstruct(T)
    error('This method does not yet support incomplete/sparse tensors.');
end

% Check the options structure.
if nargin < 2, options = struct; end
if ~isfield(options,'Range'), options.Range = 10; end
if ~isfield(options,'XMultiplier'), options.XMultiplier = 1; end

% Compute set of truncations based on upper bound on the error.
N = ndims(T);
size_tens = size(T);
maxrank = prod(size_tens)/max(size_tens);
size_tens(size_tens > maxrank) = maxrank;
in = arrayfun(@(i)1:i,size_tens,'UniformOutput',false);
out = cell(1,N);
[out{:}] = ndgrid(in{:});
out = cell2mat(cellfun(@(i)i(:),out,'UniformOutput',false));
[U,S,sv] = mlsvd(T);
loss = cellfun(@(i)[0;cumsum(flipud(i(:).^2))],sv,'UniformOutput',false);
relerr = zeros(size(out,1),1);
for n = 1:N
    lossn = loss{n};
    relerr = relerr+lossn(size_tens(n)-out(:,n)+1);
end
relerr = min(sqrt(relerr)/norm(sv{1}),1);

% Compute front of Pareto optimal truncations using a simple cull method.
mem = @(s)(sum(bsxfun(@times,size(T),s),2)+prod(s,2))/numel(T);
x = mem(out); y = relerr;
d = sqrt(x.^2+y.^2);
[d,idx] = sort(d); out = out(idx,:); x = x(idx); y = y(idx);
p = 1;
while sum(p) ~= length(x)
    idx = (x >= x(p) & y > y(p)) | (x > x(p) & y >= y(p));
    if any(idx)
        out = out(~idx,:); x = x(~idx); y = y(~idx); d = d(~idx);
    end
    p = p+1-sum(idx(1:p));
end
[~,idx] = sortrows([mem(out) -y]);
out = out(idx,:); x = x(idx); y = y(idx);

% Compute L-curve corner using the triangle method [1].
logx = options.XMultiplier*x(1:end-1); logy = y(1:end-1);
ab = cat(3,bsxfun(@minus,logx,logx.'),bsxfun(@minus,logy,logy.'));
ac = cat(3,logx(end)-logx.',logy(end)-logy.');
area = bsxfun(@times,ab(:,:,1),ac(:,:,2)) - ...
       bsxfun(@times,ab(:,:,2),ac(:,:,1));
cosa = bsxfun(@times,ab(:,:,1),ac(:,:,1)) + ...
       bsxfun(@times,ab(:,:,2),ac(:,:,2));
cosa = bsxfun(@rdivide,cosa./sqrt(ab(:,:,1).^2+ab(:,:,2).^2), ...
                       sqrt(ac(:,:,1).^2+ac(:,:,2).^2));
cosa(area >= 0 | tril(true(size(cosa))) | cosa <= cos(7*pi/8)) = -1;
[a,c] = max(max(cosa));
if a == -1, c = length(x); end

% Display output.
size_core = out(c,:);
L_ub = [out y];
if nargout == 0
    % Plot the Pareto front.
    memout = mem(out);
    style = {'Marker','+','MarkerSize',2.5};
    loglog(memout(c),y(c),'LineStyle','none','Marker','s'); hold on;
    loglog(memout,y,'b',style{:}); hold off;
    ylabel('frob(lmlrares(T,U,S))/frob(T)'); ylim([y(end-1),y(1)]);
    xlabel('(numel(U)+numel(S))/numel(T)'); xlim([memout(1) memout(end)]);
    set(datacursormode(gcf),'UpdateFcn',@datacursor);    
end
if (nargout == 0 || nargout == 3) && options.Range >= 0
    % Compute the T-MLSVD error in a range around the L-curve corner.
    idx = max(1,c-options.Range):min(length(y)-1,c+options.Range);
    relerr = zeros(size(idx));
    for i = 1:length(idx);
        Ui = cell(1,N);
        for n = 1:N, Uin = U{n}; Ui{n} =  Uin(:,1:out(idx(i),n)); end
        s = arrayfun(@(i)1:i,out(idx(i),:),'UniformOutput',false);
        relerr(i) = frob(T-lmlragen(Ui,S(s{:})))/norm(sv{1});
    end
    if nargout == 0
        hold on; loglog(mem(out(idx,:)),relerr,'r',style{:}); hold off;
        legend(['L-curve corner at ' mat2str(out(c,:))], ...
               'Upper bound on error','T-MLSVD error','Location','SW');
    end
    L_tmlsvd = [out(idx,:) relerr(:)];
end

function txt = datacursor(~,event_obj)
    pos = get(event_obj,'Position');
    txt = {['Core tensor size: ' mat2str(out(memout == pos(1),:))], ...
           ['Compression: ' num2str(100*pos(1)) '%'], ...
           ['Relative error: ' num2str(pos(2))]};
end

end
