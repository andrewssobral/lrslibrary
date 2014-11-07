function Energy = GCO_Swap(Handle,MaxIter)
% GCO_Swap   Run alpha-beta-swap algorithm.
%    GCO_Swap(Handle) runs alpha-beta-swap until convergence.
%    GCO_Swap(Handle,MaxIter) runs at most MaxIter swap cycles.
%    Returns the energy of the computed labeling.
%    The labeling itself can be retrieved via GCO_GetLabeling.
%    The order of expansion can be influenced by GCO_SetLabelOrder.
%
%    Note that neither label costs nor sparse data costs are currently
%    implement for alpha-beta-swap.

GCO_LoadLib();
if (nargin < 1), error('Swap requires handle to GCO instance'); end
if (nargin < 2), MaxIter = 1000000; end
Energy = gco_matlab('gco_swap',Handle,int32(MaxIter));
end
