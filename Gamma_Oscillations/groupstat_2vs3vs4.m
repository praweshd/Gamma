close all
clear all
clc
%%
 
pgd2 = dir('*_PGD2x2*');
cluster2 = load(pgd2.name);

 
pgd3 = dir('*_PGD.mat');
cluster3 = load(pgd3.name);

 
pgd4 = dir('*_PGD4x4*');
cluster4 = load(pgd4.name);

 
pgd5 = dir('*_PGD5x5*');
cluster5 = load(pgd5.name);
  
pgd6 = dir('*_PGD6x6*');
cluster6 = load(pgd6.name);


%% Extract percentage sig of 2x2 

for i = 1 : length(cluster2.Cluster_res) 
    PGD = cluster2.Cluster_res(i).PGD ;        
    medianval2(i) = median(PGD);       
end 

for i = 1 : length(cluster3.Cluster_res) 
    PGD = cluster3.Cluster_res(i).PGD ;        
    medianval3(i) = median(PGD);       
end 

for i = 1 : length(cluster4.Cluster_res) 
    PGD = cluster4.Cluster_res(i).PGD ;        
    medianval4(i) = median(PGD);       
end 

for i = 1 : length(cluster5.Cluster_res) 
    PGD = cluster5.Cluster_res(i).PGD ;        
    medianval5(i) = median(PGD);       
end 
for i = 1 : length(cluster6.Cluster_res) 
    PGD = cluster6.Cluster_res(i).PGD ;        
    medianval6(i) = median(PGD);       
end 


%%

close all
parulac = parula(10);

figure(1)
scatter(15*ones(length(cluster2.Cluster_res),1),medianval2,'MarkerFaceColor',parulac(1,:),'MarkerEdgeColor',parulac(1,:),'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.6)
hold on
scatter(30*ones(length(cluster3.Cluster_res) ,1),medianval3,'MarkerFaceColor',parulac(2,:),'MarkerEdgeColor',parulac(2,:),'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.6)
hold on
scatter(45*ones(length(cluster4.Cluster_res) ,1),medianval4,'MarkerFaceColor',parulac(3,:),'MarkerEdgeColor',parulac(3,:),'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.6)
hold on
scatter(60*ones(length(cluster5.Cluster_res) ,1),medianval5,'MarkerFaceColor',parulac(4,:),'MarkerEdgeColor',parulac(4,:),'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.6)
hold on
scatter(75*ones(length(cluster6.Cluster_res) ,1),medianval6,'MarkerFaceColor',parulac(5,:),'MarkerEdgeColor',parulac(5,:),'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.6)
xlabel('Clusters 2x2 -> 6x6')
ylabel('Median PGD for all clusters')
%%
threshx = 5:60;
threshy = PGD_thresh.*ones(56,1);

plot(threshx, threshy,'Color','r','LineStyle','-','Linewidth',0.4)
legend('2x2 cluster','3x3 cluster','4x4 Cluster','Shuffle Threshold')


xlim([5 60])
ylim([0 1])
ylabel('Median PGD')
print('3isbetterthan24', '-depsc', '-tiff','-r300', '-painters') 



