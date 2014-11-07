function Handle = GCO_Create(NumSites,NumLabels)
% GCO_Create    Create a GCoptimization object.
%    Handle = GCO_Create(NumSites,NumLabels) creates a new GCoptimization
%    object and returns a 'handle' to uniquely identify it.
%    Call GCO_Delete(Handle) to delete the object and free its memory.
%    Call GCO_Delete(GCO_ListHandles) to delete all GCO objects.

GCO_LoadLib();
if (nargin < 2), error('Expected 2 arguments'); end
Handle = gco_matlab('gco_create_general',int32(NumSites),int32(NumLabels));
end
