function G =  getGroupNonOverlapColor(row, col)


   % N = (row-2)*(col-2);
 N =  row * col ;
    g = sparse(zeros(N,1));
    g = diag(g);
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
    
    G = sparse(3*N, 3*groupCount);

    
    G(1:N, 1:groupCount) = g;
    G(N+1:2*N, groupCount+1:2*groupCount) =g;
    G(2*N+1:3*N, 2*groupCount+1:3*groupCount) =g;
    
end
