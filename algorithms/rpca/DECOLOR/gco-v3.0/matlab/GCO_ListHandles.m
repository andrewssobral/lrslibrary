function Handles = GCO_ListHandles()
% GCO_ListHandles     Retrieve handles to all current GCO instances
%    Useful for cleaning up GCO instances that are using memory,
%    particularly when a script was interrupted.
%    Example:
%        GCO_Delete(GCO_ListHandles);  % delete all GCO instances

GCO_LoadLib();
Handles = gco_matlab('gco_listhandles');
end
