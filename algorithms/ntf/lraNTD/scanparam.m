function varargout = scanparam(defoptions,options)
% Only proper fileds will be transferred.
% Warnings occur when options contains unexpected filed names or data type.
allfields = fieldnames(options);
opts = defoptions;
for k = 1:numel(allfields)
    if isfield(defoptions,allfields{k})&&...
            strcmp(class(options.(allfields{k})),class(defoptions.(allfields{k})))
            if ~isempty(options.(allfields{k}))
                opts.(allfields{k}) = options.(allfields{k});
            end
    else
        fprintf(strcat('Warning! Unexpected field name or data type:  [',' ',inputname(2),'.',allfields{k},'].\n'));
        fprintf('Warning! The corresponding value is not transfered.\n');
    end
end
if nargout > 1
    varargout = struct2cell(opts);
else
    varargout{1} = opts;    
end
end