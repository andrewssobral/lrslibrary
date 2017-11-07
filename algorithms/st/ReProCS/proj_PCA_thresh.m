function Q = proj_PCA_thresh(X, thresh)
    %%Function to perform the projection PCA step proposed in ReProCS paper
    %%It essentially computes the basis for the orthogonal compliment of
    %%the data matrix D on the basis matrix P of rank r.
    
    %%%             Inputs:                     %%%
    %%% X - data matrix                         %%%
    %%% thresh - threshold                      %%%
    
    
    %%%             Output:                           %%%
    %%% Q - eigenvectors with evals more than thresh  %%%
    [~, r] = size(X);
    
    [U, S, ~] = svd(X, 0);
    num_comp = min(round(r/4), length(find(diag(S) >= thresh))); 
    Q = U(:, 1 : num_comp);
end
