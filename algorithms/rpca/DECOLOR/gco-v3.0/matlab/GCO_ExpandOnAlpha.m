function GCO_ExpandOnAlpha(Handle,Alpha)
% GCO_ExpandOnAlpha   Perform a single alpha-expansion step.
%    GCO_Expansion(Handle,Alpha) takes the current labeling and performs 
%    a single expansion step on label Alpha.

GCO_LoadLib();
gco_matlab('gco_alphaexpansion',Handle,int32(Alpha));
end
