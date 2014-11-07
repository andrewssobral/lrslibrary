function Labeling = GCO_GetLabeling(Handle,varargin)
% GCO_GetLabeling     Retrieve the current labeling
%     GCO_GetLabeling(Handle) returns a column vector of all labels.
%     GCO_GetLabeling(Handle,i) returns the label of site i.
%     GCO_GetLabeling(Handle,i,count) returns labels i..i+count-1

GCO_LoadLib();
Start = int32(1);
Count = gco_matlab('gco_getnumsites',Handle);
if (length(varargin) > 2)
    error('Too many input arguments.');
end
if (length(varargin) >= 1), Start = int32(varargin{1}); Count = int32(1); end
if (length(varargin) == 2), Count = int32(varargin{2}); end
Labeling = gco_matlab('gco_getlabeling',Handle,Start,Count);
end
