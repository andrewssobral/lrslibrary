function Energy = GCO_Expansion(Handle,MaxIter)
% GCO_Expansion   Run alpha-expansion algorithm.
%    GCO_Expansion(Handle) minimizes the current energy via 
%    alpha-expansion until convergence. 
%    GCO_Expansion(Handle,MaxIter) runs at most MaxIter expansion 
%    Returns the energy of the computed labeling.
%    The labeling itself can be retrieved via GCO_GetLabeling.
%    The order of expansion can be influenced by GCO_SetLabelOrder.
%    If GCO_SetNeighbors is not called (i.e. no smoothness terms), then 
%    Expansion will internally use a greedy algorithm (no graph cuts).
%    
%    IMPORTANT: the first version uses "adaptive cycles" (changed labels) 
%    until convergence whereas the second applies up to MaxIter 
%    "standard cycles" (all labels). Each strategy is faster/slower for
%    different applications, so see what works fastest for yours.
%    

GCO_LoadLib();
if (nargin < 1), error('Expansion requires handle to GCO instance'); end
if (nargin < 2), MaxIter = -1; end
Energy = gco_matlab('gco_expansion',Handle,int32(MaxIter));
end
