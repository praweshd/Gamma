%% in paint test

function [cluster_set_interp, removerow, removecol] = inpaint3to4(mm)

rowg = 3; colg = 3; 
clusval = 0;

colnum = 9;
rownum = 13;

totclusrow = rownum-rowg+1; 
totcluscol = colnum -colg+ 1 ;

for R = 1 : totclusrow  
  for C = 1 : totcluscol
      
      clusval = clusval+1; 
      d1 = ((colnum*R-colnum)+C):((colnum*R-colnum)+C+2);
      a2 = (colnum*R+C):(colnum*R+C+2); 
      d3 = (colnum*(R+1)+C):(colnum*(R+1)+C+2); 
      
      cluster_set(clusval,:) = [d1 a2 d3]; 
      
  end      
end 

[clusrow, cluscol] = size(cluster_set);


lastcolumn_set = (colnum - 2):(colnum - 2):(rownum - 2)*(colnum - 2); %the clusters that are on the right edge
lastcolumn_set = lastcolumn_set(1:end-1); 
lastrow_set = clusrow - totcluscol + 1 : clusrow - 1;  %the clusters on the bottom edge

%If it is in the far right edge
 
if sum(ismember(mm,lastcolumn_set)) == 1

    row1 = [(cluster_set(mm,1)-1) cluster_set(mm,(1:3))];
    row2 = [(cluster_set(mm,4)-1) cluster_set(mm,(4:6))];
    row3 = [(cluster_set(mm,7)-1) cluster_set(mm,(7:9))];
    row4 = (cluster_set(mm,7)-1+cluscol) : (cluster_set(mm,7)-1+rowg+cluscol);

    %what rows and col to remove after inpaint 4x4 to get to 3x3
    removerow = 4; removecol = 1;

elseif sum(ismember(mm,lastrow_set)) == 1

    row1 = [(cluster_set(mm,1)-colnum) :(cluster_set(mm,1)-colnum+3)];
    row2 = [cluster_set(mm,(1:3)) (cluster_set(mm,3)+1)];
    row3 = [cluster_set(mm,(4:6)) (cluster_set(mm,6)+1)];
    row4 = [cluster_set(mm,(7:9)) (cluster_set(mm,9)+1)];

    removerow = 1; removecol = 4;
    
elseif mm == clusrow  

    row1 = [(cluster_set(mm,1)-colnum-1) : (cluster_set(mm,1)-colnum-1+3)];
    row2 = [ (cluster_set(mm, 1)-1)  cluster_set(mm,(1:3))];
    row3 = [(cluster_set(mm, 4 )-1)  cluster_set(mm,(4:6)) ];
    row4 = [ (cluster_set(mm, 7 )-1)  cluster_set(mm,(7:9))];

    removerow = 1; removecol = 1;

else 

    row1 = [cluster_set(mm,(1:3)) (cluster_set(mm,(3))+1)];
    row2 = [cluster_set(mm,(4:6)) (cluster_set(mm,(6))+1)];
    row3 = [cluster_set(mm,(7:9)) (cluster_set(mm,(9))+1)];
    row4 = [(cluster_set(mm,7)+colnum) :(cluster_set(mm,7)+colnum+rowg) ];

    removerow = 4; removecol = 4; 

%                          
end

 cluster_set_interp = [row1 row2 row3 row4];
 
end 

 