function [Energy D S L] = GCO_ComputeEnergy(Handle)
% GCO_ComputeEnergy   Run alpha-expansion algorithm.
%    E = GCO_ComputeEnergy(Handle) returns energy of current labeling. 
%    [E D S L] = GCO_ComputeEnergy(Handle) also provides a breakdown of
%    the energy into Data, Smooth, and Label costs.

GCO_LoadLib();
[Energy D S L] = gco_matlab('gco_computeenergy',Handle);
end
