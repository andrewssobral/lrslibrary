function G =  getGroupNonOverlapColor27(row, col)

    N =  row * col ;
    g = sparse(N,N);
 
    groupCount = 0;
   
    
    %% build non-overlapping group
    for i=2:3:col-1    
        for j=2:3:row-1  
            groupCount = groupCount+1;
            
            g(groupCount, (i-2)*row+j-1 ) = 1;
            g(groupCount, (i-2)*row+j   ) = 1;
            g(groupCount, (i-2)*row+j+1 ) = 1;
            g(groupCount, (i-1)*row+j-1 ) = 1;
            g(groupCount, (i-1)*row+j   ) = 1;
            g(groupCount, (i-1)*row+j+1 ) = 1;             
            g(groupCount, i*row+j-1 ) = 1;
            g(groupCount, i*row+j   ) = 1;
            g(groupCount, i*row+j+1 ) = 1;      

        end
    end
   
    g= g';
   
 
    g =g(:, 1:groupCount);
    
    G = [g;g;g];
    
end
