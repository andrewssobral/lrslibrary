function G = getGroupSuperColor27(I, param)
    

    N = size(I,1) * size(I, 2);

    g = sparse(N, N);
  
    imlab = vl_xyz2lab(vl_rgb2xyz(I)) ;
    imlab = single(imlab);
    slicParam = param.superpixel.slicParam;
     
    groupCount = 0;
    for i=1:size(slicParam,1)
                  
        segments = vl_slic(imlab, slicParam(i,1),  slicParam(i,2)) ;
        segments = segments(:);
       
        maxLabel = max(segments);
        for j = 0:maxLabel
            g(:, groupCount+j+1) = segments==j;
        end
        groupCount  = groupCount+ maxLabel+1;
        
        
    end
    
    % croping G
    g =g(:, 1:groupCount);
     
    G = [g;g;g];
    
end