function GCO_SetLabeling(Handle,Labeling)
% GCO_SetLabeling     Sets the current labeling
%     GCO_SetLabeling(Handle,Labeling) sets the entire labeling.

GCO_LoadLib();
if (isnumeric(Labeling))
    NumSites = gco_matlab('gco_getnumsites',Handle);
    NumLabels = gco_matlab('gco_getnumlabels',Handle);
    if (length(Labeling) ~= NumSites)
        error('Labeling must be of length NumSites');
    end
    if (~isa(Labeling,'int32'))
        if (any(floor(Labeling) ~= Labeling))
            error('Labeling was not integer valued');
        end
        Labeling = int32(Labeling);
    end
    if (min(Labeling) < 1 || max(Labeling) > NumLabels)
        error('Label must be in range 1..NumLabels');
    end
    gco_matlab('gco_setlabeling',Handle,Labeling);
end
