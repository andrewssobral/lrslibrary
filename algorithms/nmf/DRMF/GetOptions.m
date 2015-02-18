function [varargout]=GetOptions(options, varargin)
%[varargout]=GetOptions(options, varargin)
%process options given in a struct
%usage: [param_1, ...]=GetOptions(options_struct, 'param_1', param_1_default_value, ...)
% author: Liang Xiong (lxiong@cs.cmu.edu)

if isempty(options)
    options = struct('no_way_you_can_use_this',[]);
end

if nargin == 2
    opNames = varargin{1};
    if ~iscellstr(opNames); error('parameter list not right'); end

    n = length(opNames);
    if nargout > n; error('too many output'); end
    
    for ind=1:n
        varargout{ind}=options.(opNames{ind});
    end
else
    list = varargin;
    if mod(length(list), 2) ~= 0; error('parameter input not right'); end
    
    opNames = list(1:2:end);
    opDefs = list(2:2:end);
    if ~iscellstr(opNames); error('parameter list not right'); end

    n = length(opNames);
    if nargout > n; error('too many output'); end
    
    varargout = opDefs;
    for ind=1:n
        if isfield(options, opNames{ind})
            varargout{ind}=options.(opNames{ind});
        end
    end
end
