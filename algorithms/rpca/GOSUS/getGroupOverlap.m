function g =  getGroupOverlap(row, col)


    N = row*col;

    g = sparse(zeros(N,1));
    g = diag(g);
 
 
    
    %% build overlapping group
    
    % top let corner group
     
    g(1, 1 ) = 1;
    g(1, 2 ) = 1;
    g(1, 1+row ) = 1;
    g(1, 2+row ) = 1;
    
     % bottom let corner group 
  
    g(row, row-1 ) = 1;
    g(row, row ) = 1;
    g(row, row-1+row ) = 1;
    g(row, row+row ) = 1;
    
    % top right corner group
   
    g((col-1)*row+1, (col-1)*row+1 ) = 1;
    g((col-1)*row+1, (col-2)*row+1 ) = 1;
    g((col-1)*row+1, (col-1)*row+2 ) = 1;
    g((col-1)*row+1, (col-2)*row+2 ) = 1;
    
    % bottom right corner group
    
    g((col-1)*row+row , (col-1)*row+row-1 ) = 1;
    g((col-1)*row+row , (col-2)*row+row-1 ) = 1;
    g((col-1)*row+row , (col-1)*row+row ) = 1;
    g((col-1)*row+row , (col-2)*row+row ) = 1;
            
            
    
    
    % boundary group
    
    for i=2:col-1
        % top row rgoup
 
        j = 1;    
  
        g( (i-1)*row+j  , (i-2)*row+j   ) = 1;
        g( (i-1)*row+j  , (i-2)*row+j+1 ) = 1;
         
        g( (i-1)*row+j  , (i-1)*row+j   ) = 1;
        g( (i-1)*row+j  , (i-1)*row+j+1 ) = 1;             
        
        g( (i-1)*row+j  , i*row+j   ) = 1;
        g( (i-1)*row+j  , i*row+j+1 ) = 1;

        % bottom row group
    
        j = row;
        
        g( (i-1)*row+j  , (i-2)*row+j-1 ) = 1;
        g( (i-1)*row+j  , (i-2)*row+j   ) = 1;
      
        g( (i-1)*row+j  , (i-1)*row+j-1 ) = 1;
        g( (i-1)*row+j  , (i-1)*row+j   ) = 1;
                
        g( (i-1)*row+j  , i*row+j-1 ) = 1;
        g( (i-1)*row+j  , i*row+j   ) = 1;
         
        
        
    end
    
    for j=2:row-1
          % left column rgoup
    
            
        i=1;
        
        g( (i-1)*row+j  , (i-1)*row+j-1 ) = 1;
        g( (i-1)*row+j  , (i-1)*row+j   ) = 1;
        g( (i-1)*row+j  , (i-1)*row+j+1 ) = 1;             
        g( (i-1)*row+j  , i*row+j-1 ) = 1;
        g( (i-1)*row+j  , i*row+j   ) = 1;
        g( (i-1)*row+j  , i*row+j+1 ) = 1;

        % right  column group
 

        i=col;
        
        
        g( (i-1)*row+j  , (i-2)*row+j-1 ) = 1;
        g( (i-1)*row+j  , (i-2)*row+j   ) = 1;
        g( (i-1)*row+j  , (i-2)*row+j+1 ) = 1;
        g( (i-1)*row+j  , (i-1)*row+j-1 ) = 1;
        g( (i-1)*row+j  , (i-1)*row+j   ) = 1;
        g( (i-1)*row+j  , (i-1)*row+j+1 ) = 1;             
        
    end
    
    for i=2:col-1    
        for j=2:row-1  
            
            
            g( (i-1)*row+j  , (i-2)*row+j-1 ) = 1;
            g( (i-1)*row+j  , (i-2)*row+j   ) = 1;
            g( (i-1)*row+j  , (i-2)*row+j+1 ) = 1;
            g( (i-1)*row+j  , (i-1)*row+j-1 ) = 1;
            g( (i-1)*row+j  , (i-1)*row+j   ) = 1;
            g( (i-1)*row+j  , (i-1)*row+j+1 ) = 1;             
            g( (i-1)*row+j  , i*row+j-1 ) = 1;
            g( (i-1)*row+j  , i*row+j   ) = 1;
            g( (i-1)*row+j  , i*row+j+1 ) = 1;      

        end
    end
   


end