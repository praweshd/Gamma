%% PGD2_plotdata - after implementation of NAN Method 4 _ REGULAR
% This program creates two figures.
% 1. ROSE PLOT with travelling direction of the gamma
% 2. 11 x 7  grid showing where the gamma shows a travelling behavior. 
% Last modified: 9/19/2018 
%Revised July 12, 2019 for Gamma Detection PGD   
%

%%
clear all
close all 
clc
%%
%Load the states
state_mat = dir('*-states*');
load (state_mat.name);
StateIntervals = ConvertStatesVectorToIntervalSets(states);                 % 6 Intervalsets representing sleep states
REM = StateIntervals{5}; 
NREM = StateIntervals{3};
WAKE = StateIntervals{1};
state = NREM;

%Load LFP Filename
filename = dir('*.lfp');
fbasename = regexprep(filename.name,'.lfp','');

%Load Gamma Res File 
resfile = dir('*_highGAM_RES*');
load(resfile.name); 

%Load Bad channels
badch = dir('*bad_CH*');
load(badch.name);
%%
%PGDs
% allpgd = dir('*_PGD4x4*');
% load(allpgd.name);
% gridsize = 4; 
% savename_fig = [fbasename, '_ROSE4']; 

%%
%PGDs
allpgd = dir('*_PGD3x3*');
load(allpgd.name);
gridsize = 3; 
savename_fig = [fbasename, '_ROSE3']; 

%%
% %Load Shuffle Threshold mat created - THESE ARE UPLOADED IN GOOGLE DRIVE -Technical Session Dec 7 Folder
% shuff = dir('*_Shuffle*');
% load(shuff.name); 
%%

%Load the mean shuffle value for this patient's grid
%PGD_thresh = SHUFFLE(1).meanTH;
PGD_thresh = 0.1;

%  
F_in = filename.name;
CH_noise = bad_ch'; 
CH_Nall=xml2CH_N(cat(2, F_in(1:end-4),'.xml'));

Rs=1250;                                                %Sampling Frequency
duration = 0.5*Rs;
                                                        %PGD threshold based on shuffling 
yay = 0; 

%Load grid info 
NG = load('/home/prawesh/Documents/Insync/MATLAB_codes/Gamma_Oscillations/NG2019_3_info.mat');  
NGmap = NG.NG_2019_v3.map; 
NGdim = NG.NG_2019_v3.dim;

NREM_dur = length(find(states==3)); 

maxxr = 100; step = maxxr/2; 

%% 3x3 
[clusrow, cluscol, cluster_set] = makeclus(gridsize);
colg = gridsize; rowg = gridsize;
empsp = CH_noise; 
subplR = NGdim(2)-(rowg-1);
subplC = NGdim(1)-(rowg-1);
%%
close all
clc
fig=figure_ctrl('Gamma_PGD',600,1000);
 

for s = 1 : clusrow    
        
    PGD = Cluster_res(s).PGD ;  %IEDcoupled_ or IEDnotcoupled
    dirc = Cluster_res(s).angles; %in radians 

    largePGD = find(PGD>PGD_thresh);     %PGD trials greater than 0.6
    PGDsig = PGD(largePGD);              %values of PGD


    direcsig = dirc(largePGD);           %angle of direction vector for large PGDs 
        
    subaxis(subplR,subplC,s,'SpacingVert',0.01,'SpacingHoriz', 0.01, 'ML', 0.01, 'MR',0.01, 'MT', 0.01, 'MB', 0.01)
    polarhistogram(direcsig,30) 
    pax = gca;
    %rticks([0:step:maxxr]); 
    rlim([0  maxxr])
    thetaticks([0 90 180 270])
    pax.ThetaDir = 'clockwise';
    pax.ThetaZeroLocation = 'right';        
 

        if s ~= 0

            set(gca,'ThetaTickLabel',[])             

        end 
        
        if s ~=1 
            
            set(gca,'RTickLabel',[])
            
        end 

        if s == clusrow
%             rticks([0:step:maxxr])
            title([regexprep(fbasename,'_',' ')])
        end 
        
end 

%
% savename_rose = [fbasename, '_ROSE']; 
% print(savename_rose,'-dpng')

%print(savename_fig, '-depsc', '-tiff','-r300', '-painters')

print(savename_fig, '-dpng')
%print(savename_fig,'-depsc', '-tiff','-r300', '-painters')

%%
% 
% 
% 
% 
% %% Find empty spindle res
% 
% j=0;
% tempsp = [];
% for i = 1:64
%     if length(Spindles_only(i).res) == 0        
%         j = j+1;
%         tempsp(j) = i;
%         
%     end     
% end  
% 
% CH_noise = bad_ch'; 
% baddch = CH_noise(CH_noise<65);
% empsp = [ tempsp baddch']; 
% 
% 
% %% Thresholding PGD 
% allmean = cell(64,1);
% 
% binnum = 10; 
% window = 3;  
% 
% fig=figure_ctrl('PGD_Rose',1500,1000);
% title([fbasename '_ROSE'])
%  
% for s = 1 : clusrow    
%      
%         tcase = ismember(cluster_set(s,:),empsp);  
% 
% 
%         PGD = Cluster_res(s).PGD ;  %IEDcoupled_ or IEDnotcoupled
%         dirc = Cluster_res(s).angles; %in radians 
%          
% 
%         %% Extract only PGD > threshold
% 
%         
%         largePGD = find(PGD>PGD_thresh);     %PGD trials greater than 0.6
%         PGDsig = PGD(largePGD);              %values of PGD
%         
%         percent = 100*(length(largePGD)/length(PGD));
% 
%           
%         Cluster_res(s).sigper = percent; 
%        
%         
%         direcsig = dirc(largePGD);           %angle of direction vector for large PGDs 
%         
%                 
%         ang_tot = sort(rad2deg(direcsig));         
%         setN = ang_tot(ang_tot < 0); 
%         setS = ang_tot(ang_tot > 0);
%         
%         degdir_hist = hist(ang_tot,20);
% 
%         window = 3;                         %number of bins for convolution
%         alpha = 0.05;
% 
%         [dumy, pred, dumy ] = cch_conv(round(degdir_hist),window);
%         hiBound = poissinv( 1 - alpha, pred);
%         loBound = poissinv( alpha, pred);
%         hiB = hiBound;
%         loB = loBound;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%        %Plot full histogram
%         
% %         figure(66);
% %         subplot(3,1,1)    
% %         bar(degdir_hist)
% %         hold on
% %         plot(hiB,'--r')
% %         plot(loB,'--r')
% %         %ylim([0 80])
% %         
% %         subplot(3,1,2)
% %         bar(hist(setN,binnum))
% %         hold on
% %         plot(hiB(1:10),'--r')
% %         plot(loB(1:10),'--r')
% %         %ylim([0 80])
% %         
% %         subplot(3,1,3)
% %         bar(hist(setS,binnum))
% %         hold on
% %         plot(hiB(11:20),'--r')
% %         plot(loB(11:20),'--r')
% %         %ylim([0 80])
% % 
% % 
% %         figure(111)
% %         h = polarhistogram(direcsig,40);
% %         pax = gca;
% %         pax.ThetaDir = 'clockwise';
% %         pax.ThetaZeroLocation = 'right';
% %         thetaticks([0 90 180 270])
% %         rlim([0  35])
% %         rticks([0 20 35])
% %         set(gca, 'Color', 'none')
% %         transname = ['clus', num2str(s),'.png']
% %         %export_fig(transname, '-transparent')
%         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %%
% 
%         
%         subaxis(sizgrid,sizgrid,s,'SpacingVert',0.01,'SpacingHoriz', 0.01, 'ML', 0.02, 'MR',0.02, 'MT', 0.03, 'MB', 0.03)
%         polarhistogram(direcsig,30) 
%         pax = gca;
%         rticks([0:step:maxxr])
%         thetaticks([0 90 180 270])
%         pax.ThetaDir = 'clockwise';
%         pax.ThetaZeroLocation = 'right';        
%         rlim([0  maxxr])
%         
%         if s ~= 0
% 
%             set(gca,'ThetaTickLabel',[])             
% 
%         end 
%         
%         if s ~=1 
%             
%             set(gca,'RTickLabel',[])
%             
%         end 
% 
%         if s == 1
%             rticks([0:step:maxxr])
%             title([regexprep(fbasename,'_',' ')])
%         end 
%         
%      
%       
%  
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%         
%         r = 50;
%  
%         
%         for kkk = 1:2
%             
%             if kkk == 1
%                 degdir = hist(setN,binnum);
%                 test = sort(setN);
%                 hiBn = hiB(1:10);
%                 loBn = loB(1:10);
%             else
%                 degdir = hist(setS,binnum);
%                 test = sort(setS);
%                 hiBn = hiB(11:20);
%                 loBn = loB(11:20);
%             end
% 
%         
%            % Scan histogram that crosses hiB
%             for i = 1 : (length(degdir))
% 
%                 if degdir(i)>  hiBn(i)
%                     scan(i) = degdir(i);
%                 else 
%                     scan(i) = 0; 
% 
%                 end
% 
%             end
% 
%             [row,col,vbin] = find(scan');
%             
% 
%             if isempty(row) == 0
%             k=0; j=0;
%                 for i = 1 : length(row)
%                     a1 = row(i);
%                     if a1 ==1 
%                         a2 = 1;
%                     else 
%                         a2 = sum(degdir(1:a1-1));
% 
%                     end   
% 
%                     if a2 ~= 0 
%                         val = vbin(i);
%                         result = mean(test(a2:a2+val));        
%                         meanval(i) = deg2rad(result);
%                 
%                     end 
% 
%                 end 
% 
%             else
%             meanval = 0; r = 0;
% 
%             end
% 
%             mng = rad2deg(meanval);
%             [valt, idx] = max(vbin);
%             mu = mng(idx);
%             
%             
%             if length(idx) == 0 
%                 mu = 0; u1 = 0; v1 = 0; r=0;
%             end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %           
%             main(kkk) = mu; 
%             
% %             figure(2)
% %             subaxis(sizgrid,sizgrid,s)
% %             rose(direcsig,40)
% %             view([0 -90])
% %             hold on
% %             x = 0; y = 0; 
% %             r = (mu*(1/mu))*50;
% %             u1 = r * cos(deg2rad(mu)); v1 = r * sin(deg2rad(mu));
% %             quiver(x,y,u1,v1,'linewidth',2,'Color', 'red')
% % %             title([num2str(s)])
%             
% 
%         end
%         
%         allmean{s} = main;
%         
%         %Angles values in degrees 
%         Cluster_res(s).arrowang = main; 
%         
%         %Maximum Bin Value 
%         Cluster_res(s).maxbin = max(degdir_hist); 
%         
%         %XYZ Coordinates
%         xyzval = str2num(Cluster_res(s).clus);
%         Cluster_res(s).centerMNI = IED(xyzval(5)).XYZ; 
%         
% %         
% %     else 
% %     
% %     main = []; 
% %  
% %     Cluster_res(s).arrowang = main; 
%     
% %     end 
%     
%     
%     
% end 
% 
% %% PLOT AND SAVE ROSE
% % 
% %savename_rose = [fbasename, 'ROSE_overall.eps']; 
% %print(savename_rose, '-depsc', '-tiff','-r300', '-painters')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %% PLOT AND SAVE BIN HEAT MAP
% close all
% 
% %Load all IED channels
% %% 
% % allied = dir('*PGD_SPIED_COUPLING*');
% % load(allied.name);
% % 
% % %% find IED
% % 
% % [r,c]=size(PGD_IED_coupling);
% % for i=1:r
% %    IED_CHs(i)=PGD_IED_coupling(i,1).IED_ch;
% % end
% 
% %%
% close all
% 
% for i = 1 : clusrow
%     
%         %tcase = ismember(cluster_set(i,:),empsp);
%     
% %     if sum(tcase) ~= 0  
% %         
% %         percentsig(i) = NaN; 
% %         
% %     else
%         
%         binval(i) = Cluster_res(i).sigper;
%         
% %     end 
% 
% end
% 
% savename_bin = [fbasename , '_PGDbin_overall.eps']; 
% 
% % cmin = (min(binval)+(min(binval)/2));
% % cmax =max(binval)-50;
% 
% cmin = min(binval)+5;
% cmax =max(binval)-5;
% 
% fig=figure_ctrl('Sig%',700,700);
% binval = binval';
% sigrid = reshape(binval,6,6)';
% imagesc(sigrid); colormap jet;  colorbar; caxis ([cmin cmax ])
% title([fbasename,' PGD bin heat'],'Interpreter', 'none')
% 
% dq = interp2(sigrid,8);
% 
% fig=figure_ctrl('Sig%',700,700);
% imagesc(dq); colormap jet;  colorbar; caxis ([cmin cmax])
% title(['interp8'])
% 
% 
% %down to 8x8 
% for i = 1:8 
%     
%     %colp = 160*(i-1)+1;
%     colp = 160*i;
%     new1grid(:,i) = dq(:,colp);
%     
% end 
% 
% for i = 1 :8
%     
%     %rowp = 160*(i-1)+1;
%     rowp = 160*i;
%     newgrid(i,:) = new1grid(rowp,:);
% end 
% % 
% fig=figure_ctrl('Fin%',700,700);
% imagesc(newgrid); colormap jet;  colorbar; caxis ([cmin cmax])
% title('Final')
% savename_binheat = [fbasename , '_PGD_sigperheat.eps']; 
% print(savename_binheat, '-depsc', '-tiff','-r300', '-painters') 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %print(savename_bin, '-depsc', '-tiff','-r300', '-painters')
% 
% 
% %% SAVE THE FINAL MAT FILE
% 
% 
% savename = [fbasename, '_PGD_NAN_M4.mat']; 
% save(savename,'Cluster_res')
% 
% 
% 
% 
% 
