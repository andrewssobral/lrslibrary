function GCO_Delete(Handle)
% GCO_Delete    Delete a GCoptimization object.
%    GCO_Delete(Handle) deletes the object corresponding to Handle 
%    and frees its memory.

gco_matlab('gco_delete',int32(Handle));
end
