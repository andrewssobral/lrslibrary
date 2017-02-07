function [res loc] = minmaxk(mexfun, list, k, dim, varargin)
% function [res loc] = minmaxk(mexfun, list, k, dim)
% 
% Return in RES the K smallest/largest elements of LIST
%   RES is sorted in ascending/descending order
% [res loc] = minmaxk(...)
%   Location of the smallest/largest: RES=LIST(LOC)
% [res loc] = minmaxk(..., dim)
%   specify the dimension to operate
% [res loc] = minmaxk(..., dim, 'sorting', false)
%   to disable the post-sorting step (true by default)
%
% Author Bruno Luong <brunoluong@yahoo.com>
% Contributor: Matt Fig (suggest return same loc as SORT for full blowed
%                        result)
% Last update: 24/May/2009: work on sparse matrix list
%              28/Jun/2009: used inplace columns for full-array
%              10/Aug/2009: releaseinplace is called to cleanup
%                           if MEX fails
%              10/Jan/2010: possibility to disable post-sorting step

clist=class(list);
% Mex functions requires input in double
if ~strcmpi(clist,'double')
    list=double(list);
end

% Look for single selection value by default
if nargin<3 || isempty(k)
    k=1;
else
    k=double(k);
end

szlist = size(list);
if nargin<4
    if isvector(list) && szlist(1)==1
        dim=2;
    else
        dim=1;
    end
end

if mod(length(varargin),2)
    error('MINMAXK: options must come as property/value pairs');
end

postsorting = getoptionpair({'postsorting', 'sorting', 'sort'}, ...
                            true, varargin);

nd=ndims(list);
if dim<1 || dim>nd
    error('MINMAXK: dim must be between 1 and ndims(LIST)=%d', nd);
end

% Will be used for sorting
if isequal(mexfun,@minkmex)
    smode='ascend';
else
    smode='descend';
end

% Do we need to get location? 
getloc=nargout>=2;

% Put operating dimension to the first
list=shiftdim(list,dim-1);

% operating length
l=size(list,1);
% Number of vectors
szl=size(list);
N=prod(szl(2:end));

szres=szl;
k=min(k,l);
szres(1)=k;
res=zeros(szres,clist); % Allocate array having the same class with list
if getloc
    loc=zeros(szres,'double');
end
if k>=l % Matt Fig's suggestion
    res = list; 
    if getloc
        repvec=size(loc); repvec(1)=1;
        loc = repmat((1:k)',repvec);
    end
else
    try % use try/catch instead of onCleanup for better compatibility
        if getloc
            if ~verLessThan('MATLAB', '8.3') || issparse(list) 
                for n=1:N
                    [res(:,n) loc(:,n)] = mexfun(list(:,n),k); %#ok
                end
            else
                for n=1:N
                    cn = inplacecolumnmex(list,n); % inplace column
                    [res(:,n) loc(:,n)] = mexfun(cn,k);
                    releaseinplace(cn);
                    %[res(:,n) loc(:,n)] = mexfun(list(:,n),k);
                end
            end
        else
            if ~verLessThan('MATLAB', '8.3') || issparse(list) 
                for n=1:N
                    res(:,n) = mexfun(list(:,n),k);
                end
            else
                for n=1:N
                    cn = inplacecolumnmex(list,n); % inplace column
                    res(:,n) = mexfun(cn,k);
                    releaseinplace(cn);
                    %res(:,n) = mexfun(list(:,n),k);
                end
            end
        end
    catch
        % If something is wrong
        % It crashes if cn is not released properly
        if exist('cn','var') && ~isempty(cn)
            releaseinplace(cn);
        end
        % rethrow the error (likely memory)
        rethrow(lasterror);
    end
end

% This is the post processing step of sorting the selection data
% The purpose is to have a nicely formatted output, that's all

if getloc
    if postsorting
        [res is] = sort(res,1,smode);
        j=(0:N-1)*k;
        % Use reshape instead of bsxfun for backward compatible
        if exist('bsxfun','builtin')
            is = bsxfun(@plus, reshape(is,[k N]), j);
        else
            is = reshape(is,[k N]) + repmat(j,[k 1]);
        end
        loc = reshape(loc(is),size(loc));
    end
    % Put the operating dimension at the right place
    loc = shiftdim(loc,nd+1-dim);
else
    if postsorting
        res = sort(res,1,smode);
    end
end

% Put the operating dimension at the right place
res = shiftdim(res,nd+1-dim);

end % minmaxk

function val = getoptionpair(name, defaultval, vargin)
% Get the value from property/value pairs
    val = defaultval;
    for k=1:2:length(vargin)
        if strmatch(vargin{k},name)
            val = vargin{k+1};
            return
        end
    end
end % getoptionpair