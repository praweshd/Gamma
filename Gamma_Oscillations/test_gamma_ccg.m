%% Testing CCG between detected Gamma events
% Prawesh Dahal
%June 24

close all
clear all
clc
%%
good_CH = [0,1,2,3,4,6,20,21,23,25,27,29,31,33,34,35,36,37,38,39,41,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114];
%% createbadchannels([])
createbadchannels([setdiff(1:128,good_CH+1)]-1)
%%
%Load the states
state_mat = dir('*-states*');
load (state_mat.name);
StateIntervals = ConvertStatesVectorToIntervalSets(states);                 % 6 Intervalsets representing sleep states
REM = StateIntervals{5}; 
NREM = StateIntervals{3};
WAKE = StateIntervals{1};
state = NREM;
pp = pwd; 
ani = pp(43:49);
%%
%Load CH
lfp_file = dir('*.lfp');
[~, fbasename, ~] = fileparts(lfp_file.name);
CH_Nall=xml2CH_N(cat(2,lfp_file.name(1:end-4),'.xml'));

%%
%Load detected event 
detec_file = dir('*_highGAM_RES.mat');
%detec_file = dir('*_highGAM_res_II.mat');
load(detec_file.name); 

%Load Bad channels
badch = dir('*bad_CH*');
load(badch.name);
CH_noise = bad_ch'; 

NREM_dur = length(find(states==3)); 
Rs = 1250;

%% New Grid?
%Load grid info 
NG = load('/home/prawesh/Documents/Insync/MATLAB_codes/Gamma_Oscillations/NG2019_4_info.mat'); 
 
NGmap = NG.NG_2019_v4.map; 
NGdim = NG.NG_2019_v4.dim;
%% Old Grid
%Load grid info 
NG = load('/home/prawesh/Documents/Insync/MATLAB_codes/Gamma_Oscillations/NG2019_3_info.mat'); 
 
NGmap = NG.NG_2019_v3.map; 
NGdim = NG.NG_2019_v3.dim;

%%  OCCURENCE RATE GRID PLOT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Occurence Rate  
close all
clc
for i = 1 : length(gammas) 
    
    tmp = length(gammas(i).res(:,1));     
    occ = tmp/NREM_dur;  
    if ismember(i,CH_noise) == 1      
        gammas(i).occ = NaN;   
    else 
        gammas(i).occ = occ;
    end
end  

k=0; 
for i = 1:NGdim(2)    
    for j = 1 : NGdim(1)        
        k = k+1; 
        indx = NGmap(k); 
        gam_occ(i,j) = gammas(indx).occ;       
        
    end     
end 

fig=figure_ctrl('Occ Rate',400,600);
imagesc(gam_occ); 
colormap jet; colorbar; 
% caxis([0.5 1.3])
title(['Occurence Rate']) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% savename_fig = [fbasename, '_H_occrate.eps']; 
% print([fbasename, '_H_occrate'], '-depsc', '-tiff','-r300', '-painters')

print([fbasename, '_H_occrate'], '-dpng')

%% -----TRIGGER AVERAGE-------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------
%% Trigger Average Power
close all
clc
Rs = 1250;
duration = 0.2*Rs;
trial_N=500;

for i = 1  %: length(NGmap)
    
    if ismember(NGmap(i),CH_noise)
        disp(['No trigger average for bad ch ', num2str(NGmap(i))])
        i = i +1; 
    else 
        close all
        disp(['Running Trigger Average on CH ', num2str(NGmap(i))])

        gam = gammas(NGmap(i)).res;

        if length(gam) < 500
            trial_N = length(gam); 
        end 

        %time from res OR
        avg_res = (gam(:,1) + gam(:,3))/2;
        time = round(avg_res*Rs);
        time = time(1:trial_N);
         

        %time from NScope
        % time=round(3929.470*Rs);

        data = Dat_tracker(lfp_file.name,time,duration,CH_Nall);
        data(CH_noise,:,:) = NaN;
        data=data(NGmap,:,:); CH_N = length(NGmap); 

        %Reshape 
        tmp=reshape(data,CH_N,trial_N*duration);

        %Filter Hilbert
        n=3; Wn=[60 80];
        [b,a]=butter(n,2*Wn/Rs,'bandpass');
        Data_fil_tmp=filtfilt (b,a,tmp')';
        Data_hil_tmp=hilbert(Data_fil_tmp')';

        %Reshape back
        data_fil=reshape(Data_fil_tmp,CH_N,duration,trial_N);
        data_hil_all=reshape(Data_hil_tmp,CH_N,duration,trial_N);  

        data_hil = abs(data_hil_all);

        midphiT = round(length(data(1,:,1))/2);

        %data_hil_phase = unwrap(angle(data_hil_all(:,midphiT,:)));
        %Trigger Average
        data_hil_avg=mean(data_hil,3);

        %Time window for max power
        trigger_window = 50; %how many ms around the time of the midpoint of event 
        trigger_window = (trigger_window/1000)*Rs;

        start_point=round((duration/2) - trigger_window );
        end_point= round(start_point   + 2*trigger_window);

        %Power
        gam_hil=data_hil_avg(:,start_point:end_point);
        gam_pwr = max(gam_hil,[],2); 
        %gam_pwr=sum(gam_hil,2)./(length (start_point:end_point)); 
        %gam_pwr=sum(gam_hil,2)./median(gam_hil,2);

        grid_pwr=reshape(gam_pwr,9,13);
%%
%         close all

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Grid Filter and Hilbert
        fig=figure_ctrl('power',900,1300);
        for j=1:117
            subplot(13,9,j);
            trial = 1; 
            %plot (data(j,start_point:end_point, trial),'k'); hold on 
            plot (data_fil(j,start_point:end_point, trial),'r'); hold on 
            plot (data_hil(j,start_point:end_point, trial)+50,'k'); hold on 
            ylim([-100 200])
            axis_cleaner

            if j== find(NGmap == i)
               title(['CH', num2str(i), ' T= ', num2str(trial), 'res=', num2str(time(trial))]) 
            end 

        end
%         print([fbasename, 'trigavgstep'], '-depsc', '-tiff','-r300', '-painters')

        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fig=figure_ctrl('grid power',600,800);
        imagesc(grid_pwr'); colormap jet;  
        colorbar; 
        caxis([50 100])
        title(['High Gamma Trig Avg CH', num2str(NGmap(i))]) 
        % 

        
%        print([fbasename, '_tgravg', num2str(NGmap(i))],'-dpng')        
%        print([fbasename, 'trigavgstepFIG'], '-depsc', '-tiff','-r300', '-painters')


    end

end



%% -----CCG------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------

%% CCG per trial

close all
clc 
fig=figure_ctrl('CCG_BAR',1000,2000);

for ref = 68 %:  length(NGmap)
    
    refch = NGmap(ref);
    res1= ((gammas(refch).res(:,1)) + (gammas(refch).res(:,3)))./2;
    disp(['CCG for CH ', num2str(refch)])

    for chs = 1:length(NGmap)
        
        ref2ch = NGmap(chs); 
        res2= ((gammas(ref2ch).res(:,1)) + (gammas(ref2ch).res(:,3)))./2;

        window=500;
        bin=10;

        [H, B, hiB, loB]=CCG_res(res1,res2,window,bin,Rs);

        subaxis(13,9,chs, 'SpacingVert',0.01,'SpacingHoriz', 0.01, 'ML', 0.01, 'MR',0.01, 'MT', 0.01, 'MB', 0.01)
        bar(B,H); 
        % title('HC-Cx ripple CCG')
        hold on; plot(B,hiB,'--r','LineWidth',0.5); 
        hold on; plot(B,loB,'--r','LineWidth',0.5); 
%         % win=5000;
        ylim([0 30]);
        xlim([-400 400])
%         % title (num2str(RIPPLE_DET.ECx(CH).CH))
          
%         GammaCCG(ref).zerobin(chs) =  H(B==0);

    end 

end 

print(['CCGbar_', num2str(refch)], '-depsc', '-tiff','-r300', '-painters') 

%% Regular CCG
close all
clc 
for ref = 1 :  length(gammas)

    res1= ((gammas(ref).res(:,1)) + (gammas(ref).res(:,3)))./2;
    disp(['CCG for CH ', num2str(ref)])

    for chs = 1:length(gammas)
        
        res2= ((gammas(chs).res(:,1)) + (gammas(chs).res(:,3)))./2;

        window=500;
        bin=10;

        [H, B, hiB, loB]=CCG_res(res1,res2,window,bin,Rs); 
        GammaCCG(ref).zerobin(chs) =  H(B==0);

    end
end 

%%
gamma_file = strcat(fbasename, '_CCG'); 
save(gamma_file, 'GammaCCG');

%% Plot CCG, create global CCG
 
ccgfile = dir('*_CCG.mat');
load(ccgfile.name);

%CH_noise = []; bad_ch =[];

close all
clc
for i = 1 : length(GammaCCG)
    if ismember(i,CH_noise) || ~ismember(i,NGmap)
        disp(['No CCG for bad ch ', num2str(i)])
        i = i +1;                              
    else 
        
        close all
        ccgval = GammaCCG(i).zerobin';
        ccgval(bad_ch) = NaN; 
        ccgval = ccgval(NGmap,1);                                                    
        grid_ccgX=reshape(ccgval,NGdim(1),NGdim(2));        
        
        fig=figure_ctrl('grid ccg',600,800);
        heatmap(grid_ccgX','CellLabelColor','none');   
%         imagesc(grid_ccgX')
        colormap jet;  colorbar; 
        caxis([min(ccgval)  max(ccgval)])
        title(['High GammaCCG ', num2str(i)])
         
        print([fbasename, '_CCG_', num2str(i)], '-dpng')
%         print(['CCGheat_', num2str(i)], '-depsc', '-tiff','-r300', '-painters') 


    end
end

 %%
%Global
close all
clc
ccgfile = dir('*_CCG.mat');
load(ccgfile.name);

globaltest = zeros(length(GammaCCG),1);

for i = 1:length(GammaCCG) 
    indccg = GammaCCG(i).zerobin';
    indccg(i) = 0;
    globaltest = globaltest + indccg;   
end 

globaltest(bad_ch) = NaN; 
gridtestt = globaltest(NGmap,1);
grid_gob = reshape(gridtestt,NGdim(1),NGdim(2));

gridtestnorm = normalize(gridtestt,'range');
grid_gobnorm = reshape(gridtestnorm,NGdim(1),NGdim(2));

ccgnorm.val = [gridtestnorm]; 
ccgnorm.sum = nansum(gridtestnorm)
ccgnorm_file = strcat(fbasename, '_CCG_VAL'); 
save(ccgnorm_file, 'ccgnorm');

fig=figure_ctrl('glob ccg',600,800); 
heatmap(grid_gob','CellLabelColor','none');  
colormap jet;  colorbar; 
caxis([0 5000])
title([ani,' ',fbasename, ' Global'])
print([ani,' ',fbasename, '_Global'], '-dpng')

fig=figure_ctrl('glob ccg norm',600,800); 
heatmap(grid_gobnorm','CellLabelColor','none');  
colormap jet;  colorbar; 
caxis([0 1])
title([ani,' ',fbasename, ' Global'])
print([ani,' ',fbasename, '_GlobalNORM'], '-dpng')
 
%% -----PGD CCGs-----------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------
%% PGD Parameters
%% Parameters for PGD
%PGDs
allpgd = dir('*_PGD3x3*');
load(allpgd.name);

% %Load Shuffle Threshold mat created - THESE ARE UPLOADED IN GOOGLE DRIVE -Technical Session Dec 7 Folder
% shuff = dir('*_Shuffle*');
% load(shuff.name); 
% PGD_thresh = SHUFFLE(1).meanTH;
PGD_thresh = 0.4; 
 
%%Create a 3x3 cluster list of electrodes 
gridsize = 3; 
[clusrow, cluscol, cluster_set] = makeclus(gridsize);
colg = gridsize; rowg = gridsize;
empsp = CH_noise; 


%% Check cross-correlation between two clusters that have good modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clc
load('/home/prawesh/Documents/Insync/MATLAB_codes/Gamma_Oscillations/3x3subplot.mat')

left = find(yoyo==0);
right = find(yoyo==1);
 
avgdelay = 0;

for cA = 1:77   
    close all
    figure_ctrl('grid ccg travel',2000,2000);
    clusresA = Cluster_res(cA).RES;
    
    PGDA = Cluster_res(cA).PGD;   
    clusresA = clusresA(PGDA>PGD_thresh); 
    
    for cB = 1:77

        clusresB = Cluster_res(cB).RES;
        PGDB = Cluster_res(cB).PGD;
        clusresB = clusresB(PGDB>PGD_thresh); 
        window=700;
        bin=10;

        [H, B, hiB, loB]=CCG_res(clusresA,clusresB,window,bin,Rs);
        
         subplot(11,11,left(cB))

         if cA==cB
            maxind = find(H==max(H));            
            b = bar(B,H);
            b.FaceColor = 'flat';
            b.CData(maxind(1),:) = [0.6350 0.0780 0.1840];

        else 
            H(find(B==0)) = 0; 
            b = bar(B,H);             
            maxind = find(H==max(H)); 
            b.FaceColor = 'flat';
            b.CData(maxind(1),:) = [0.9290 0.6940 0.1250];
         
        end 
        hold on; plot(B,hiB,'--r'); 
        hold on; plot(B,loB,'--r');
        %set(gca,'YTickLabel',[]);
        xlim([- 200 200]);
        ylim([0 5])
        
        %NET DELAY CALCULATION
        zeroindx = find(B==0);
        posms = find(B==100);
        prems = find(B==-100);
        Hpos = sum(H(zeroindx+1 : posms));
        Hpre = sum(H(prems : zeroindx-1));        
        netdelay(cB) = Hpos - Hpre;         

    end
     
     avgdelay = avgdelay + normalize(netdelay,'range'); 

     subplot(11,11,right)
     netdelaygrid = reshape(netdelay',7,11)';         
     imagesc(netdelaygrid); colormap jet;  colorbar; caxis([-max(netdelay) max(netdelay)]) 
     title(['Delay map clus ', num2str(cA)])
     print([fbasename, '_CCGdelay_', num2str(cA)], '-dpng')    
end



close all
fig=figure_ctrl('grid ccg pgd',600,800);
avgdelayN = avgdelay./77; 
avgdelaygrid = reshape(avgdelayN',7,11)';         
imagesc(avgdelaygrid); colormap jet;  colorbar; caxis([0 .8])
title('AVERAGE DELAY FROM 77 TRIALS NORMALIZED')
print([fbasename, '_AVGdelay_'], '-dpng') %%
% print([fbasename, '_AVERAGEdelay_ALLCLUS', num2str(77)],'-depsc', '-tiff','-r300', '-painters')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% CCG DELAY IED- Jose Jennifer
close all
clear all
clc
%%

HPch = 167; 
mpfc = [130 131 133 134];
iedstr = dir('*_K16_post*');
savename1= ['per_NRAT_C19_K16_HP_', num2str(HPch)];
savename2= ['del_NRAT_C19_K16_HP_', num2str(HPch)];

load(iedstr.name);

close all
clc
load('plotmap.mat') 

left = find(plotmap==0);
right = find(plotmap==1);
mpfcloc = find(plotmap==2); 
NG = load('/home/prawesh/Documents/Insync/MATLAB_codes/Gamma_Oscillations/NG2019_3_info.mat'); 
 
NGmap = NG.NG_2019_v3.map; 
NGdim = NG.NG_2019_v3.dim;

allCH = [IED_all.channel];
HPres = IED_all(find(allCH == HPch)).res;
Rs = 1250;
figure_ctrl('grid ccg travel',2000,2000);

for i = 1 : length(NGmap)
    
    REFCH = NGmap(i);
    whereCH = find(allCH == REFCH);
    
    subplot(13,15,left(i))
    if isempty(whereCH)
        B=0;
        H=0;
        b = bar(B,H);
        netdelay(i) = NaN; sigper(i) = NaN; sigval(i) = NaN;    
    
    elseif ~isempty(whereCH)
        NGres = IED_all(whereCH).res;
        window=500;
        bin=10;        
        if isempty(NGres)
            B=0;
            H=0;
            b = bar(B,H);
            netdelay(i) = NaN; sigper(i) = NaN; sigval(i) = NaN;            

        else
            [H, B, hiB, loB]=CCG_res(HPres,NGres,window,bin,Rs);
            b = bar(B,H);
            hold on; plot(B,hiB,'--r'); 
            hold on; plot(B,loB,'--r');
            %set(gca,'YTickLabel',[]);
            xlim([- 100 100]); ylim([0 5])

%             %NET DELAY CALCULATION
%             zeroindx = find(B==0);
%             posms = find(B==30);
%             prems = find(B==-30);
%             Hpos = sum(H(zeroindx+1 : posms));
%             Hpre = sum(H(prems : zeroindx-1));        
%             netdelay(i) = Hpos - Hpre;            
            
             %Zero-bin percentage
            cchEvt = (H.*length(HPres)*bin/1000)';

            if hiB(find(B==0)) == 0
                sigper(i) = 0 ;
            elseif hiB(find(B==0)) ~= 0
                if H(find(B==0)) > hiB(find(B==0))

                    sigper(i) = cchEvt(find(B==0))/length(NGres);
                else
                    sigper(i) = 0 ;
                end 
            end
          
            %New net delay
            tmp = intersect(find(B>0), find(B<30));
            if sum(hiB) ~=0
                tmp = intersect(tmp, find(H>hiB));
                 
                sigval(i) = sum(H(tmp));  
            elseif sum(hiB) == 0
                sigval(i) = 0; 
            end
        end  
    end 
end 

sigpergridm = reshape(sigper',9,13)';
sigvalgridm = reshape(sigval',9,13)';


subplot(13,15,right)
heatmap(sigvalgridm); colormap jet;  colorbar; %caxis([-2  2]) 
% netdelaygrid = reshape(netdelay',9,13)';         
% heatmap(netdelaygrid,'CellLabelColor','none'); colormap jet;  colorbar; caxis([-2 2])% caxis([-max(netdelay) max(netdelay)]) 
% title([savename]) 
 
for i = 1:length(mpfc)    
    
    whereCH = find(allCH == mpfc(i));
    mpfcres = IED_all(whereCH).res;
    
    subplot(13,15,mpfcloc(i))
    [H, B, hiB, loB]=CCG_res(HPres,mpfcres,window,bin,Rs);
    b = bar(B,H);
    hold on; plot(B,hiB,'--r'); 
    hold on; plot(B,loB,'--r');
    %set(gca,'YTickLabel',[]);
    xlim([- 100 100]);
    ylim([0 5])

    %NET DELAY CALCULATION
    zeroindx = find(B==0);
    posms = find(B==100);
    prems = find(B==-100);
    Hpos = sum(H(zeroindx+1 : posms));
    Hpre = sum(H(prems : zeroindx-1));        
    netdelaym(i) = Hpos - Hpre;
    
    %Zero-bin percentage
    cchEvt = (H.*length(HPres)*bin/1000)';
    if H(find(B==0)) > hiB(find(B==0))
        sigperm(i) = cchEvt(find(B==0))/length(mpfcres);
    else 
        sigperm(i) = 0; 
    end

    %New net delay
    tmp = intersect(find(B>0), find(B<30));
    tmp = intersect(tmp, find(H>hiB)); 
    
    sigvalm(i) = sum(H(tmp));             
    
end
 
subplot(13,15,mpfcloc(5:8))
% netdelaygridm = reshape(netdelaym',4,1)';
sigvalgridm = reshape(sigvalm',4,1)';

% heatmap(netdelaygridm','CellLabelColor','none'); colormap jet;  colorbar; caxis([-2  2]) 
heatmap(sigvalgridm','CellLabelColor','none'); colormap jet;  colorbar; %caxis([-2  2]) 
print([savename2], '-dpng') 
 
figure_ctrl('grid ccg',1000,1000);
heatmap(sigpergridm); colormap jet;  colorbar; caxis([0 0.3])  
print([savename1], '-dpng') 
title('No. of sig zero bin IED events/total IED events')

%%
%close all

REFCH = NGmap(33);
whereCH = find(allCH == REFCH);
NGres = IED_all(whereCH).res;
[H, B, hiB, loB]=CCG_res(HPres,NGres,window,bin,Rs);
b = bar(B,H);
hold on; plot(B,hiB,'--r'); 
hold on; plot(B,loB,'--r');
cchEvt = (H.*length(HPres)*bin/1000)';
%caxis([-max(netdelay) max(netdelay)]) 
%print([savename], '-dpng') 

 
%% -----CCGs of two direction gamma -----------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------
%% Any different in spatial extent of one of other direction travelling gamma



%%
close all
clc
maxxr = 300;
step = maxxr/2; 
for cA = 1:77

    clusresA = Cluster_res(cA).RES;
    dircA = Cluster_res(cA).angles; %in radians 
    PGDA = Cluster_res(cA).PGD;   

    direcsigA = dircA(PGDA>PGD_thresh); direcsigA = sort(direcsigA); 
    clusresA = clusresA(PGDA>PGD_thresh); clusresA = sort(clusresA); 


    dirX = direcsigA( 1 : (round(length(direcsigA)/2)) ); 
    dirY = direcsigA( (round(length(direcsigA)/2)+1) : end ); 

    resX = clusresA( 1 : (round(length(direcsigA)/2)) ); 
    resY = clusresA( (round(length(direcsigA)/2)+1) : end  ); 

    res1 = resX; 

    for chs = 1:128

        res2= ((gammas(chs).res(:,1)) + (gammas(chs).res(:,3)))./2;
        window=500;
        bin=10;
        [H, B, hiB, loB]=CCG_res(res1,res2,window,bin,Rs);
        zerobinX(chs) =  H(B==0);

    end 
 
    ccgval = zerobinX';
    ccgval = ccgval(NGmap,1);
    grid_ccgX=reshape(ccgval,9,13);

    res1 = resY; 

    for chs = 1:128

        res2= ((gammas(chs).res(:,1)) + (gammas(chs).res(:,3)))./2;
        window=500;
        bin=10;
        [H, B, hiB, loB]=CCG_res(res1,res2,window,bin,Rs);
        zerobinY(chs) =  H(B==0);

    end 

    ccgval = zerobinY';
    ccgval = ccgval(NGmap,1);
    grid_ccgY=reshape(ccgval,9,13);

    close all
    fig=figure_ctrl('grid ccg travel',600,800);

    subaxis(2,2,1,'SpacingVert',0.01,'SpacingHoriz', 0.01, 'ML', 0.02, 'MR',0.02, 'MT', 0.03, 'MB', 0.03)
    polarhistogram(dirX,30) 
    pax = gca;
    rticks([0:step:maxxr]); 
    rlim([0  maxxr])
    thetaticks([0 90 180 270])
    pax.ThetaDir = 'clockwise';
    pax.ThetaZeroLocation = 'right';  
    title(['Cluster ' num2str(cA)])

    subaxis(2,2,2,'SpacingVert',0.01,'SpacingHoriz', 0.01, 'ML', 0.02, 'MR',0.02, 'MT', 0.03, 'MB', 0.03)
    polarhistogram(dirY,30) 
    pax = gca;
    rticks([0:step:maxxr]); 
    rlim([0  maxxr])
    thetaticks([0 90 180 270])
    pax.ThetaDir = 'clockwise';
    pax.ThetaZeroLocation = 'right';    

    subaxis(2,2,3,'SpacingVert',0.01,'SpacingHoriz', 0.01, 'ML', 0.02, 'MR',0.02, 'MT', 0.03, 'MB', 0.03)
    imagesc(grid_ccgX'); colormap jet;  colorbar; caxis([0 35])

    subaxis(2,2,4,'SpacingVert',0.01,'SpacingHoriz', 0.01, 'ML', 0.02, 'MR',0.02, 'MT', 0.03, 'MB', 0.03)
    imagesc(grid_ccgY'); colormap jet;  colorbar; caxis([0 35])
    
    savename_fig = [fbasename, '_pgdccg', num2str(cA)];
    print(savename_fig, '-dpng')

end
 