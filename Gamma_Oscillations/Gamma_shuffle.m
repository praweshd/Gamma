%% Shuffle and find 95% threshold for each cluster. 
%Be in the folder of the animal/patient
% Prawesh Dahal
%Created November 30, 2018
%Finalized December 14, 2018
%Revised July 11, 2019 for Gamma Detection PGD   


%% LOAD THE WAVELET SPINDLE RES FILE



%%
close all
clear all
clc

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
%%
%Load Bad channels
badch = dir('*bad_CH*');
load(badch.name);

% %Load Shuffle Threshold mat created - THESE ARE UPLOADED IN GOOGLE DRIVE -Technical Session Dec 7 Folder
% shuff = dir('*_Shuffle*');
% load(shuff.name); 
% 
% %Load the mean shuffle value for this patient's grid
% PGD_thresh = SHUFFLE(1).meanTH;
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

%% Create a 3x3 cluster list of electrodes 

rowg = 3; colg = 3; 
j = 0;

colnum = 9;
rownum = 13;

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
empsp = CH_noise; 
 

%% Find available workable clusters 
  
for mm = 1 : clusrow
    
   CH = cluster_set(mm,:);      
    if sum(ismember(empsp,NGmap(CH))) == 0         
        salvageclus(mm) = mm;        
    else
        salvageclus(mm) = NaN; 
    end    
end 

salvageclus = reshape(salvageclus, 7, 11)'; 

%% Manually choosing five cluster to shuffle  

shufclus = [7 40 49 57 77];


%%

for i = 1 : length(shufclus)  
    cls = shufclus(i); 
    newcluster(i,:) = cluster_set(cls,:);
end 


%% Find threshold for those five clusters 

[newclusrow, newcluscol] = size(newcluster);

for mm = 1: newclusrow
    
    disp(['PPU working on cluster shuffle ', num2str(mm)]) 
 
    clear samplernd
    clear spi_sort
    clear SPI_time
    clear trialS
    trialS = []; 
    spi_sort = 0;
    SPI_time = 0;
    samplernd = cell(1);
    
    for i = 1:newcluscol    
        CH = NGmap(newcluster(mm,i));

        if ismember(empsp,CH) ~= 1

            gam = gammas(CH).res;
            %[row, col] = size(spindles); 
            %avg_spt = gam(:,1); 
            avg_spt = (gam(:,1) + gam(:,3))/2;
            timernd{i} = avg_spt;
            samplernd{i} = round(timernd{i}*Rs);

        end

    end
    
    spi_sampc = reshape(samplernd,[],1);
    spi_samp = cell2mat(spi_sampc); 
    spi_sort = sort(spi_samp);

%%

     if length(spi_sort) > 0 
         
        disp('in the loop')    
        SPI_time = []; j=0; 
        
        for i = 1 : (length(spi_sort)-1)

            A = spi_sort(i+1); 
            B = spi_sort(i);
            tmp = spi_sort(i+1) - spi_sort(i);

            if tmp >= 0.2*Rs
                j = j+1;
                SPI_time(j) = spi_sort(i);
            end       

        end 

        SPI_time = sort(SPI_time'); 
        
     else
         SPI_time = [];
        
     end
    
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVES CONCATENATED SPINDLES TIMES IN THE CELL 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if length(SPI_time) > 0
        disp('SPI time big 0') 
        
        if length(SPI_time) > 800            
            disp('SPI time > 800')            
            loopn = ceil(length(SPI_time)/800);
            for i = 1: loopn 

                chk = 800*(i-1);
                tmp = length(SPI_time) -chk;

                if tmp>= 800
                    trialS(i,1) = 1+(800*(i-1));
                    trialS(i,2) = 800*i;

                elseif tmp < 800
                    trialS(i,1) = 1+(800*(i-1));
                    trialS(i,2) = 800*(i-1) + tmp; 
                end  
            end 
            
        elseif length(SPI_time) <= 800
            disp('SPI time < 800')
            trialS(1,1) = 1;
            trialS(1,2) = length(SPI_time);
        end 

        clear trial_start
        clear trial_end
        
        for qqq = 1 : length(trialS(:,1))      
        
            trial_start = trialS(qqq,1);
            trial_end   = trialS(qqq,2); 
            
            data= Dat_tracker(F_in,SPI_time(trial_start:trial_end),duration,CH_Nall);        
            data(CH_noise,:,:) = 0;
            
            data=data(NGmap,:,:); CH_N = length(NGmap);              
            trial_N = (trial_end - trial_start + 1); % THIS IS THE TRIAL SET OF 800 

            %Reshape 
            tmp=reshape (data,CH_N,trial_N*duration);

            %Filter Hilbert
            n=3; Wn=[60 80];
            [b,a]=butter(n,2*Wn/Rs,'bandpass');
            Data_fil_tmp=filtfilt (b,a,tmp')';
            Data_hil_tmp=hilbert (Data_fil_tmp')';

            %Resample to extract phase data    
            resample_factor=1 ;
            Data_hil_tmp_resample=resample (Data_hil_tmp', 1,resample_factor)';

            data_fil=reshape(Data_fil_tmp,CH_N,duration,trial_N);
            data_hil_all=reshape(Data_hil_tmp,CH_N,duration/resample_factor,trial_N);

%             t=(1:duration/resample_factor)./(Rs/resample_factor)-0.5;
%             timezero = round(length(t+1)/2);
            midphiT = round(length(data(1,:,1))/2);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        %%
            indiv_clus = newcluster(mm,:);
            clust_hil = data_hil_all(indiv_clus, :, :); %Now extracting the hilbert data for those cluster channels 
   
            for trial = 1 : trial_N

                for ii = 1:length(indiv_clus)
                    instphase(ii,1) = unwrap((angle(clust_hil(ii,midphiT,trial))));   %Filling my good 3x3 cluster with phase values
                end
                
                gradgrid = reshape(instphase,[colg,rowg])';               
                instphase_FIN = reshape(gradgrid',[1,cluscol])';

                [px, py] = gradient(gradgrid);

                ppx = reshape(px', length(indiv_clus), 1);
                meanppx = mean(ppx);

                ppy = reshape(py', length(indiv_clus), 1);      
                meanppy = mean(ppy);

                vecmag(trial,1) = norm([meanppx meanppy]); 
                num = vecmag(trial,1); 

                for k = 1:length(indiv_clus)

                    newn(k,1) = norm([ppx(k) ppy(k)]);

                end 

                den = mean(newn);

                PGD(trial,1) = num/den; 
                dirc(trial,1) = atan2(-meanppy,-meanppx);  
                
                %%% NOW SHUFFLE 
                
                scount = 500; 
                shuf_array = instphase_FIN;
                
                    for kk = 1:scount
            
                        shuf_array = shuf_array(randperm(length(shuf_array)));
                        shufgrid = reshape(instphase_FIN(randperm(length(instphase_FIN))), [colg,rowg])';   

                        [sx, sy] = gradient(shufgrid);       
                        spx = reshape(sx', cluscol, 1);
                        spy = reshape(sy', cluscol, 1);
                        meanspx = mean(spx); meanspy = mean(spy);         
                        num_shuf = norm([meanspx meanspy]);         

                            for i = 1:cluscol

                                new_s(i,1) = norm([spx(i) spy(i)]);

                            end 

                        den_shuf = mean(new_s);                
                        PGD_shuf(kk) = num_shuf/den_shuf; 
                       
                    end
                
                    PGD_SHUF(trial,1) = median(PGD_shuf); 
                

            end 

            

            SHUFFLE(mm).clus = num2str(NGmap(newcluster(mm,:)));
            SHUFFLE(mm).ShufMedPGD((trial_start:trial_end),1) = PGD_SHUF; 
            SHUFFLE(mm).ActPGD((trial_start:trial_end),1) = PGD;
                
            clear PGD
            clear PGD_shuf
            clear PGD_SHUF
               
            
            %name1 = ['NY632_Day6_Part2_2x2_' num2str(indiv_clus(1)) 'PGD.mat']; 
            %save(name1,'PGD_fin')
            
            yay = 1;
        end 
            
    elseif length(SPI_time) <= 0 
        
        disp('SPI_time is less 0')        

            SHUFFLE(mm).clus = num2str(NGmap(newcluster(mm,:)));
            SHUFFLE(mm).ShufMedPGD((trial_start:trial_end),1) = []; 
            SHUFFLE(mm).ActPGD((trial_start:trial_end),1) = [];
                
    end
       
    
    close all

    shufflemedian = SHUFFLE(mm).ShufMedPGD;
    actual = SHUFFLE(mm).ActPGD;

    actPGD_fin = median(actual);    
    SHUFFLE(mm).medianact =  actPGD_fin;
    
    thresh = prctile(shufflemedian,95);
    SHUFFLE(mm).threshold95 =  thresh;
    
    %%%PLOT%%%%
    
%     y = 0:900;      
%     x1 = ones(length(y),1)*max(thresh); 
% 
%     figure(1)
%     fig=figure_ctrl('SHUFF',500,500);
%     histogram(shufflemedian,25,'FaceColor','b')
%     hold on
%     plot(ones(length(y),1).*actPGD_fin, y, 'Linewidth',2,'Color','k')
%     plot(x1, y, 'Linewidth',2,'Color','r')
%     xlabel('Median PGD')
%     ylabel('Counts of shuffles')
%     legend('Shuffled PGD','Actual data', '95%')
%     name = ['Shuffle', num2str(cluster_set(mm,1))];
%     %saveas(fig,name,'jpg')
    
    clear actPGD_fin
    clear thresh       
    clear x1
    
       

end
    


%% Calculate mean threshold for five cluster 

for i = 1:newclusrow 
    
    allthresh(i) = SHUFFLE(i).threshold95;   
    
    
end 

meanTH = mean(allthresh); 
devTH  = std(allthresh); 

SHUFFLE(1).meanTH = meanTH; 
SHUFFLE(1).devTH  = devTH; 
SHUFFLE(1).salvageclus = salvageclus;  
SHUFFLE(1).shufclusrows = shufclus; 
    
 
 
%  
%  
%  
%% SAVE MAT FILE 


 savename = [fbasename, '_Shuffle.mat']; 
 save(savename,'SHUFFLE')
%  
% %  
% %  
% %  
% % %  %%
% % %  
% % %  
% % % %Plot 
% % % close all
% % % 
% % % shufflemedian = SHUFFLE(2).ShufMedPGD;
% % % actual = SHUFFLE(2).ActPGD;
% % % 
% % % actPGD_fin = median(actual); 
% % % y = 0:900; 
% % % thresh = prctile(shufflemedian,95);
% % % x1 = ones(length(y),1)*max(thresh); 
% % % 
% % % figure(1)
% % % histogram(shufflemedian,25,'FaceColor','b')
% % % hold on
% % % plot(ones(length(y),1).*actPGD_fin, y, 'Linewidth',2,'Color','k')
% % % plot(x1, y, 'Linewidth',2,'Color','r')
% % % xlabel('Median PGD')
% % % ylabel('Counts of shuffles')
% % % legend('Shuffled PGD','Actual data', '95%')
% % % %xlim([ 0 0.8])
% % % 
% % 
% % % for i = 1 : 27
% % %     
% % %     TEST(i) = SHUFFLE(i).medianact;
% % %     
% % %     
% % % end 
% % %      