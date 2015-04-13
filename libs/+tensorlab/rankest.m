function [R,L_lb,L_cpd] = rankest(T,options)
%RANKEST Estimate rank.
%   rankest(T) plots an L-curve of number of rank-one terms in a canonical
%   polyadic decomposition. The x-axis corresponds to the number of
%   rank-one terms, and the y-axis corresponds to the relative error of the
%   CPD in that many rank-one terms. Additionally, the corner R of the
%   resulting L-curve is estimated.
%
%   R = rankest(T) does not plot anything and instead returns the number of
%   rank-one terms R corresponding to the corner of the L-curve.
%
%   [R,L_lb,L_cpd] = rankest(T) also returns the L-curve ranks L_lb(:,1)
%   and L_cpd(:,1) and the corresponding lower bound on the truncation
%   error L_lb(:,2) and relative error of the CPD approximation L_cpd(:,2),
%   respectively.
%
%   rankest(T,options) may be used to set the following options:
%
%      options.MaxR =           - Maximum number of rank-one terms to try.
%      numel(T)/max(size_tens)
%      options.MinR = 1         - Minimum number of rank-one terms to try.
%      options.MinRelErr = 1e-2 - Determines an upper threshold for the
%                                 number of rank-one terms R to try.
%      options.Solver = @cpd    - The solver used to compute the CPD for
%                                 each number of rank-one terms R. Called
%                                 as options.Solver(T,R, ...
%                                 options.SolverOptions), where
%                                 options.SolverOptions is an options
%                                 structure passed to the solver.
%      options.XMultiplier = 1  - The importance of the number of rank-one
%                                 terms in determining the L-curve corner,
%                                 relative to the importance of the
%                                 relative error of the approximation.
%
%   See also mlrankest.

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
if isstruct(T), size_tens = T.size;
else size_tens = size(T); end

% Check the options structure.
if nargin < 2, options = struct; end
if ~isfield(options,'MaxR')
    options.MaxR = prod(size_tens)/max(size_tens);
end
if ~isfield(options,'MinR'), options.MinR = 1; end
if ~isfield(options,'MinRelErr'), options.MinRelErr = 1e-2; end
if ~isfield(options,'Solver'), options.Solver = @cpd; end
if ~isfield(options,'SolverOptions'), options.SolverOptions = struct; end
if ~isfield(options.SolverOptions,'Compression')
    options.SolverOptions.Compression = false;
end
if ~isfield(options.SolverOptions,'Display')
    options.SolverOptions.Display = false;
end
if ~isfield(options,'XMultiplier'), options.XMultiplier = 1; end

% Compute lower bound on truncation error.
frobT = frob(T);
if ~isstruct(T)
    [~,~,sv] = mlsvd(T);
    Rmax = max(cellfun(@length,sv))-1;
    lb = max(cell2mat(cellfun( ...
             @(s)[sqrt(flipud(cumsum(flipud(s(:).^2)))).'/frobT ...
             zeros(1,Rmax+1-length(s))],sv(:),'UniformOutput',false)));
    lb = lb(2:end);
else
    Rmax = 0;
    lb = [];
end

% Determine the range of rank-one terms to test.
if ~isempty(lb)
    R = find(lb(options.MinR:end)<=options.MinRelErr,1,'first')+ ...
        options.MinR-1;
    if isempty(R), R = max(options.MinR,Rmax+1); end
    R = R:options.MaxR;
else
    R = 1:options.MaxR;
end

% Compute the relative error for each number of rank-one terms R(r).
relerr = zeros(1,length(R));
for r = 1:length(R)
    
    % Compute the CPD in R rank-one terms.
    [U,output] = options.Solver(T,R(r),options.SolverOptions);
    if isfield(output,'Refinement') && ...
       isfield(output.Refinement,'fval')
        relerr(r) = sqrt(2*output.Refinement.fval(end))/frobT;
    elseif isfield(output,'Algorithm') && ...
           isfield(output.Algorithm,'fval')
        relerr(r) = sqrt(2*output.Algorithm.fval(end))/frobT;
    elseif isfield(output,'fval')
        relerr(r) = sqrt(2*output.fval(end))/frobT;
    else
        relerr(r) = frob(cpdres(T,U))/frobT;
    end
    
    % Stop is options.MinRelErr has been reached.
    if relerr(r) <= options.MinRelErr
        R = R(1:r);
        relerr = relerr(1:r);
        break;
    end
    
    % Update the plot.
    if r > 1 && r < length(R) && nargout == 0, Lcurve(true); end
    
end

% Compute L-curve corner using the triangle method [1].
logx = options.XMultiplier*R(:);
logy = log10(relerr(:));
if R(1) ~= 1 && ~isempty(lb)
    logx = [options.XMultiplier;logx];
    logy = [log10(lb(1));logy];
end
ab = cat(3,bsxfun(@minus,logx,logx.'),bsxfun(@minus,logy,logy.'));
ac = cat(3,logx(end)-logx.',logy(end)-logy.');
area = bsxfun(@times,ab(:,:,1),ac(:,:,2)) - ...
       bsxfun(@times,ab(:,:,2),ac(:,:,1));
cosa = bsxfun(@times,ab(:,:,1),ac(:,:,1)) + ...
       bsxfun(@times,ab(:,:,2),ac(:,:,2));
cosa = bsxfun(@rdivide,cosa./sqrt(ab(:,:,1).^2+ab(:,:,2).^2), ...
                       sqrt(ac(:,:,1).^2+ac(:,:,2).^2));
cosa(area >= 0 | tril(true(size(cosa))) | cosa <= cos(7*pi/8)) = -1;
[a,opt] = max(max(cosa(:,2:end)));
if isempty(a), opt = 1; end
if a == -1, opt = size(cosa,2)-1; end

% Display output.
L_lb = [(1:Rmax).' lb(:)];
L_cpd = [R(:) relerr(:)];
if nargout == 0, Lcurve(false); end
R = R(opt);

function Lcurve(update)
    if R(r) <= 1, return; end
    style = {'Marker','+','MarkerSize',2.5};
    if ~update
        semilogy(R(opt),relerr(opt),'rs','LineStyle','none'); hold on;
    end
    semilogy(R(1:length(relerr)),relerr,'r',style{:}); hold on;
    semilogy(1:length(lb),lb,style{:});
    semilogy([1 R(r)],options.MinRelErr*[1 1],'k:');
    text(1,options.MinRelErr,'MinRelErr','VerticalAlignment','Top');
    hold off;
    ylim([min(options.MinRelErr/(10^0.5),min(relerr(1:r))) lb(1)]);
    xlim([1 R(r)]);
    ylabel('frob(cpdres(T,U))/frob(T)'); xlabel('R');
    xt = get(gca,'XTick'); set(gca,'XTick',xt(mod(xt,1) == 0));
    if update
        if ~isempty(lb)
            legend(['CPD error (trying R = ' int2str(R(r+1)) '...)'], ...
                'Lower bound on error','Location','NE');
        else
            legend(['CPD error (trying R = ' int2str(R(r+1)) '...)'], ...
                'Location','NE');
        end
    else
        if ~isempty(lb)
            legend(['L-curve corner at R = ' int2str(R(opt))], ...
                'CPD error','Lower bound on error','Location','NE');
        else
            legend(['L-curve corner at R = ' int2str(R(opt))], ...
                'CPD error','Location','NE');
        end
    end
    set(datacursormode(gcf),'UpdateFcn',@datacursor);
    drawnow;
end

function txt = datacursor(~,event_obj)
    pos = get(event_obj,'Position');
    txt = {['R: ' int2str(pos(1))], ...
           ['relative error: ' num2str(pos(2))]};
end

end
