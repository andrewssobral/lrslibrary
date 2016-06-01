function GCO_SetLabelCost(Handle,LabelCost,LabelSubset)
% GCO_SetLabelCost   Set costs associated with using labels.
%    GCO_SetLabelCost(Handle,LabelCost) with scalar LabelCost gives the
%    same cost to all labels.
%    GCO_SetLabelCost(Handle,LabelCost) with 1xNumLabels LabelCost
%    associates cost LabelCost(k) to label k.
%    GCO_SetLabelCost(Handle,LabelCost,LabelSubset) sets the cost for using
%    at least one label mentioned in LabelSubset (i.e. LabelSubset is a 
%    vector containing label indices). The cost is paid once.
%    SetLabelCost can be called before or after Expansion.

GCO_LoadLib();
if (nargin < 2)
    error('Expected at least 2 arguments');
end
NumLabels = gco_matlab('gco_getnumlabels',Handle);
if (length(LabelCost) ~= 1 && length(LabelCost) ~= NumLabels)
    error('LabelCost must be scalar or of length NumLabels');
end
EnergyTermClass = gco_matlab('gco_get_energyterm_class');
LabelCostClass = class(LabelCost);
if (~strcmp(LabelCostClass,EnergyTermClass))
    OldLabelCost = LabelCost;
    LabelCost = cast(OldLabelCost,EnergyTermClass);
    if (length(LabelCost) > 50 || any(any(cast(LabelCost, LabelCostClass) ~= OldLabelCost)))
        warning('GCO:type',['LabelCost converted to ' EnergyTermClass]);
    end
end
if (nargin < 3) 
    gco_matlab('gco_setlabelcost',Handle,LabelCost);
else
    if (length(LabelCost) ~= 1), error('LabelCost must be scalar'); end
    if (any(LabelSubset < 1) || any(LabelSubset > NumLabels))
        error('LabelSubset must contain indices from 1..NumLabels');
    end
    if (~isa(LabelSubset,'int32'))
        if (any(any(floor(LabelSubset) ~= LabelSubset)))
            error('LabelSubset must contain integers from 1..NumLabels');
        end
        LabelSubset = int32(LabelSubset);
    end
	gco_matlab('gco_setlabelcost',Handle,LabelCost,LabelSubset);
end
end
