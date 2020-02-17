%% in paint nans for various clusters

function [cluster_set_interp, removerow, removecol]  = inpaintany(mm, gridsize)

[clusrow, cluscol, cluster_set] = makeclus(gridsize);

rowg = gridsize; colg = gridsize; 

colnum = 9;
rownum = 13;

totclusrow = rownum-rowg+1; 
totcluscol = colnum -colg+ 1 ;

gridsl = gridsize - 1;

lastcolumn_set = (colnum - gridsl):(colnum - gridsl):(rownum - gridsl)*(colnum - gridsl); %the clusters that are on the right edge
lastcolumn_set = lastcolumn_set(1:end-1); 
lastrow_set = clusrow - totcluscol + 1 : clusrow - 1;  %the clusters on the bottom edge

cluster_set_interp = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(ismember(mm,lastcolumn_set)) == 1

    for i = 1 : gridsize+1
        
        if i ~= (gridsize+1)
            row = [(cluster_set(mm,(1+(i-1)*gridsize))-1) cluster_set(mm,(1+(i-1)*gridsize) : (1+(i-1)*gridsize) + gridsize-1)];
            
        elseif i == (gridsize+1)
            
            row = [(cluster_set(mm,(1+(i-2)*gridsize))-1+colnum) : (cluster_set(mm,(1+(i-2)*gridsize))-1+colnum + gridsize)];
            
        end
        cluster_set_interp = [cluster_set_interp row];
    end
    
removerow = gridsize+1; removecol = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

elseif sum(ismember(mm,lastrow_set)) == 1    
k=0;

for i = 1 : gridsize
    
    if k==0
        row = [(cluster_set(mm,1)-colnum):(cluster_set(mm,1)-colnum + gridsize)];
        cluster_set_interp = [cluster_set_interp row];
        k=1;
    end 
    
   row = [cluster_set(mm,(1+(i-1)*gridsize) : (1+(i-1)*gridsize) + gridsize-1) (cluster_set(mm,(1+(i-1)*gridsize + gridsize-1))+1) ];
   cluster_set_interp = [cluster_set_interp row];  
        
end 
removerow = 1 ; removecol = gridsize+1;     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif mm == clusrow  
    
    k=0;

    for i = 1 : gridsize

        if k==0
            row = [(cluster_set(mm,1)-colnum-1):(cluster_set(mm,1)-colnum-1 + gridsize)];
            cluster_set_interp = [cluster_set_interp row];
            k=1;
        end 

       row = [(cluster_set(mm,(1+(i-1)*gridsize))-1) cluster_set(mm,(1+(i-1)*gridsize) : (1+(i-1)*gridsize) + gridsize-1)];
       cluster_set_interp = [cluster_set_interp row];  

    end 
 
    removerow = 1; removecol = 1;

else 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    for i = 1 : gridsize+1
        
        if i ~= (gridsize+1)
            row = [cluster_set(mm,(1+(i-1)*gridsize) : (1+(i-1)*gridsize) + gridsize-1) (cluster_set(mm,(1+(i-1)*gridsize + gridsize-1))+1)];
            
        elseif i == (gridsize+1)
            
            row = [(cluster_set(mm,(1+(i-2)*gridsize))+colnum) : (cluster_set(mm,(1+(i-2)*gridsize))+colnum + gridsize)];
            
        end
        cluster_set_interp = [cluster_set_interp row];
    end 

    removerow = gridsize+1; removecol = gridsize+1; 
                          
end
end 
