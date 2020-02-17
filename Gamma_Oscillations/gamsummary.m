%% Perform gamma analysis loading the grand struct file 
% Prawesh Dahal
% Feb 13, 2020

close all 
clear all
clc

sumfile = dir('*_summary.mat');
load(sumfile.name)

% %Old Grid
% NG = load('/home/prawesh/Documents/Insync/MATLAB_codes/Gamma_Oscillations/NG2019_3_info.mat'); 
% NGmap = NG.NG_2019_v3.map; 
% NGdim = NG.NG_2019_v3.dim;

%New Grid OR15
NG = load('/home/prawesh/Documents/Insync/MATLAB_codes/Gamma_Oscillations/NG2019_4_info.mat'); 
NGmap = NG.NG_2019_v4.map; 
NGdim = NG.NG_2019_v4.dim;

%%
close all
clc

for ses = 1: length(ORbehavior) 
    
    for chs = 1:length(ORbehavior(ses).gam)
        ccgmap = ORbehavior(ses).gam(chs).ccgheat; 

        ccgmapN = normalize(ccgmap,'range');

 
%         fig=figure_ctrl('grid ccg',600,800);
%         heatmap(grid_ccgX','CellLabelColor','none');   
%         colormap jet;  colorbar;
%         caxis([min(ccgmapN)  max(ccgmapN)])
        
        ccgsort = sort(ccgmapN, 'descend');
        if sum(isnan(ccgsort)) > 0
            tmp = sum(isnan(ccgsort));
            ccgsort = ccgsort(1+tmp:end);
        end 
        
%         figure(1)
%         plot(1:length(ccgsort),ccgsort)
%         
        auc = trapz(1:length(ccgsort),ccgsort);
        ORbehavior(ses).gam(chs).auc = auc;

    end 

end 

%%

%Total sessions for OR13 is 33
clc
area = zeros(128,33);
for i = 1:33
    
    sesnum = find([ORbehavior.sesorder] == i);
    temp = [ORbehavior(sesnum).gam.auc]';
    
    if length(temp) ~= 128
        
        temp = [temp ;zeros(128-length(temp),1)];
        area(:,i) = temp; 
        
    else 
        area(:,i) = temp; 
               
    end 
    
end

%% Plot area across days
close all
clc

mapgrid = reshape(NGmap, 9, 13)';
regionA = mapgrid(1:5, 1:5);    regions(1).num = area(regionA(:),:)'; 
regionB = mapgrid(1:5, 6:9);    regions(2).num = area(regionB(:),:)';  
regionC = mapgrid(6:9, 1:5);    regions(3).num = area(regionC(:),:)'; 
regionD = mapgrid(6:9, 6:9);    regions(4).num = area(regionD(:),:)'; 
regionE = mapgrid(10:13, 1:5);  regions(5).num = area(regionE(:),:)'; 
regionF = mapgrid(10:13, 6:9);  regions(6).num = area(regionF(:),:)';

%OR13
% bas = 6; cbd = 15; obj = 27;

%OR6
bas = 4; cbd = 16; obj = 26;

%%
%c=parula(128);
close all
clc

fig=figure_ctrl('Occ Rate',900,1000);

for i = 1:6
    
subplot(3,2,i)
test = [regions(i).num];
shadedErrorBar(1:33,mean(test,2),std(test'),'b'); 
hold on
plot(ones(70,1)*bas, 1:70)
plot(ones(70,1)*cbd, 1:70)
plot(ones(70,1)*obj, 1:70)
xlim([1 33]); ylim([1 70])
xticks([1:5:33])
xlabel('Sessions')
ylabel('Area')

end 
title('Gamma Summary OR6')
print(['OR6 Channel-Location GAM Spatial'],'-dpng')
print(['OR6 Channel-Location GAM Spatial'], '-depsc', '-tiff','-r300', '-painters') 
%%
 


