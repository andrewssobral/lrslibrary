function G = getGroupSuperNonOverlapColor27(I, param)
    

    N = size(I,1) * size(I, 2);

    g = sparse(zeros(N,1));
    g = diag(g);
 
    imlab = vl_xyz2lab(vl_rgb2xyz(I)) ;
    imlab = single(imlab);
    slicParam = param.superpixel.slicParam;
     
    groupCount = 0;
             
    segments = vl_slic(imlab, slicParam(1,1),  slicParam(1,2)) ;
    segments = segments(:);

    for j = 0:segments(N)
        g(:, groupCount+j+1) = segments==j;
    end
    groupCount  = groupCount+ segments(N)+1;

    % croping G
    g =g(:, 1:groupCount);
     
    G = [g;g;g];   
    
     
end