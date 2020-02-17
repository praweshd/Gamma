%% Give the size of the square grid (3 if 3x3, for instance), it makes the cluster

function [clusrow, cluscol, cluster_set] = makeclus(gridsize)

colnum = 9;
rownum = 13;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if gridsize == 2
        
    rowg = 2; colg = 2; 
    j = 0;

    for R = 1 : (rownum-rowg+1)     
      for C = 1 : (colnum - colg+ 1)
          j = j+1; 
          d1 = ((colnum*R-colnum)+C):((colnum*R-colnum)+C+1);
          a2 = (colnum*R+C):(colnum*R+C+1); 
          cluster_set(j,:) = [d1 a2];
      end      
    end 
    [clusrow, cluscol] = size(cluster_set);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if gridsize == 3
        
    rowg = 3; colg = 3; 
    j = 0;

    for R = 1 : (rownum-rowg+1)     
      for C = 1 : (colnum - colg+ 1) 
          j = j+1; 
          d1 = ((colnum*R-colnum)+C):((colnum*R-colnum)+C+2);
          a2 = (colnum*R+C):(colnum*R+C+2); 
          d3 = (colnum*(R+1)+C):(colnum*(R+1)+C+2); 
          cluster_set(j,:) = [d1 a2 d3]; 
      end      
    end 
    [clusrow, cluscol] = size(cluster_set);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if gridsize == 4
        
    rowg = 4; colg = 4; 
    j = 0;

    for R = 1 : (rownum-rowg+1)     
      for C = 1 : (colnum - colg+ 1) 
          j = j+1; 
          d1 = ((colnum*R-colnum)+C):((colnum*R-colnum)+C+3);
          a2 = (colnum*R+C):(colnum*R+C+3); 
          d3 = (colnum*(R+1)+C):(colnum*(R+1)+C+3);
          d4 = (colnum*(R+2)+C):(colnum*(R+2)+C+3);
          cluster_set(j,:) = [d1 a2 d3 d4]; 
      end      
    end 
    [clusrow, cluscol] = size(cluster_set);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if gridsize == 5
        
    rowg = 5; colg = 5; 
    j = 0;

    for R = 1 : (rownum-rowg+1)     
      for C = 1 : (colnum - colg+ 1) 
          j = j+1; 
          d1 = ((colnum*R-colnum)+C):((colnum*R-colnum)+C+4);
          a2 = (colnum*R+C):(colnum*R+C+4); 
          d3 = (colnum*(R+1)+C):(colnum*(R+1)+C+4);
          d4 = (colnum*(R+2)+C):(colnum*(R+2)+C+4);
          d5 = (colnum*(R+3)+C):(colnum*(R+3)+C+4);
          cluster_set(j,:) = [d1 a2 d3 d4 d5]; 
      end      
    end 
    [clusrow, cluscol] = size(cluster_set);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if gridsize == 6
        
    rowg = 6; colg = 6; 
    j = 0;

    for R = 1 : (rownum-rowg+1)     
      for C = 1 : (colnum - colg+ 1) 
          j = j+1; 
          d1 = ((colnum*R-colnum)+C):((colnum*R-colnum)+C+5);
          a2 = (colnum*R+C):(colnum*R+C+5); 
          d3 = (colnum*(R+1)+C):(colnum*(R+1)+C+5);
          d4 = (colnum*(R+2)+C):(colnum*(R+2)+C+5);
          d5 = (colnum*(R+3)+C):(colnum*(R+3)+C+5);
          d6 = (colnum*(R+4)+C):(colnum*(R+4)+C+5);
          cluster_set(j,:) = [d1 a2 d3 d4 d5 d6]; 
      end      
    end 

    [clusrow, cluscol] = size(cluster_set);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

