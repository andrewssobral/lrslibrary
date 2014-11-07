function GCO_SetNeighbors(Handle,Weights)
% GCO_SetNeighbors   Set (weighted) pairwise connectivity of all sites.
%     GCO_SetNeighbors(Handle,Weights) determines which sites are neighbors
%     and thereby have a SmoothCost associated with them. Weights is a
%     sparse NumSites-by-NumSites matrix, where Weights(i,j) > 0 indicates 
%     that sites i and j are neighbors. If Weights is a 0-1 matrix, smooth
%     costs are spatially invariant. See SetSmoothCost for more.
%
%     SetNeighbors cannot be called after Expansion. 
%     Note: only the upper-triangular area of Weights is consulted 
%           because the connectivity is undirected. 

GCO_LoadLib();
NumSites = gco_matlab('gco_getnumsites',Handle);
if (size(Weights) ~= [ NumSites NumSites ])
    error('Neighbors must be of size [ NumSites NumSites ]');
end
if (~issparse(Weights))
    if (NumSites > 100)
        warning('Sparsifying the Neighbors matrix (performance warning)');
    end
    if (~isa(Weights,'double'))
        error('Neighbors matrix must be of type double, but with integral values');
    end
    Weights = sparse(Weights);
end
gco_matlab('gco_setneighbors',Handle,Weights);
end
