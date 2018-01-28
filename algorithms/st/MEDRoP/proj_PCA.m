function Q = proj_PCA(D, P, r)
    %%Function to perform the projection PCA step proposed in ReProCS paper
    %%It essentially computes the basis for the orthogonal compliment of
    %%the data matrix D on the basis matrix P of rank r.
    
    %%%             Inputs:                     %%%
    %%% D - data matrix                         %%%
    %%% P - basis matrix                        %%%
    %%% r - rank of the output                  %%%
    
    %%%             Output:                     %%%
    %%% Q - basis for the othogonal compliment  %%%
    
    alpha = size(D, 2);

    if(~isempty(P))
        D_proj = D - (P * (P' * D));        
    else
        D_proj = D;
    end
    
    [Q, ~] = svds(1 / alpha * (D_proj * D_proj'), r);
    
end
