function GCO_SetDataCost(Handle,DataCost,Label)
% GCO_SetDataCost   Set the data cost of individual sites.
%    GCO_SetDataCost(Handle,DataCost) accepts a NumLabels-by-NumSites 
%    int32 matrix where DataCost(k,i) is the cost of assigning
%    label k to site i. In this case, the MATLAB matrix is pointed to 
%    by the C++ code, and so an int32 DataCost array is not copied.
%
%    GCO_SetDataCost(Handle,DataCost,Label) accepts a 2-by-N int32 matrix
%    of (site,cost) pairs, i.e. DataCost(2,i) is the cost of assigning 
%    Label to site DataCost(1,i). The site ids must be sorted in increasing 
%    order. All ommitted site ids are assumed to be infeasible for Label.
%    This 'sparse' version of SetDataCost allows Expansion to run much 
%    faster when labels are only feasible for a small subset of sites.
%    It is possible to assign infeasible labelings via GCO_SetLabeling,
%    but GCO_ComputeEnergy will add huge constants to represent each
%    infeasible assignment.
%
%    SetDataCost can be called repeatedly, even after Expansion. 
%

GCO_LoadLib();
if (nargin < 2), error('Expected at least 2 arguments'); end
if (~isnumeric(DataCost)), error('DataCost must be numeric'); end
if (~isreal(DataCost)), error('DataCost cannot be complex'); end
NumLabels = gco_matlab('gco_getnumlabels',Handle);
NumSites  = gco_matlab('gco_getnumsites', Handle);
if (nargin == 2)
    Label = 0; % no specific label
    if (size(DataCost) ~= [ NumLabels NumSites ])
        error('DataCost size must be [ NumLabels NumSites ]');
    end
else
    if (Label < 1 || Label > NumLabels)
        error('Label must be in range 1..NumLabels');
    end
    if (size(DataCost,1) ~= 2)
        error('Sparse DataCost must contain two rows');
    end
end
if (~isa(DataCost,'int32'))
    if (NumSites*NumLabels > 200 || any(any(floor(DataCost) ~= DataCost)))
        % warning('GCO:int32','DataCost converted to int32');
    end
    DataCost = int32(DataCost);
end
gco_matlab('gco_setdatacost',Handle,DataCost,int32(Label));
end
