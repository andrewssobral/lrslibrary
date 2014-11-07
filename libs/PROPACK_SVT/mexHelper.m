function mexHelper(varargin)
n = length(varargin);
indx = [];
for i = 1:n
    if ~isempty( varargin{i} )
        indx = [indx,i];
    end
end

mex( varargin{indx} )