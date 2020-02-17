%% Create the MAIN PGD MAT FILE WITH NAN INCLUDED USING METHOD 4
%Be in the folder of the animal/patient
% Prawesh Dahal
%Created November 30, 2018
%Finalized December 14, 2018
%Revised July 11, 2019 for Gamma Detection PGD   

close all
clear all
clc
 
%%

% {Parameters}
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

% %Load Shuffle Threshold mat
%shuff = dir('*_Shuffle*');
%load(shuff.name); 
%PGD_thresh = SHUFFLE(1).meanTH;
PGD_thresh = 0.4;

F_in = filename.name;
CH_noise = bad_ch'; 
CH_Nall=xml2CH_N(cat(2, F_in(1:end-4),'.xml'));

Rs=1250;                                                %Sampling Frequency
duration = 0.6*Rs;
 
%Load grid info 
NG = load('/home/prawesh/Documents/Insync/MATLAB_codes/Gamma_Oscillations/NG2019_3_info.mat');  
NGmap = NG.NG_2019_v3.map; 
NGdim = NG.NG_2019_v3.dim;

NREM_dur = length(find(states==3)); 

%% Create a 3x3 cluster list of electrodes 
savename = [fbasename, '_PGD3x3.mat']; 
gridsize = 3; 
[clusrow, cluscol, cluster_set] = makeclus(gridsize);
colg = gridsize; rowg = gridsize;
empsp = CH_noise; 
%% Main Loop 

for mm = 1 :clusrow               
    disp(['PPU working on cluster ', num2str(mm)]) 
    clear samplernd spi_sort SPI_time trialS
 
    trialS = [];  spi_sort = 0; SPI_time = 0;
    samplernd = cell(1);
    
    for i = 1:cluscol  
        CH = NGmap(cluster_set(mm,i));
        if sum(ismember(empsp,CH)) == 0 
            gam = gammas(CH).res;
            avg_spt = (gam(:,1) + gam(:,3))./2;
            timernd{i} = avg_spt;
            samplernd{i} = round(timernd{i}*Rs);
        end
    end 

%  Concatenate all the spindle times
    spi_sampc = reshape(samplernd,[],1); 
    spi_samp = cell2mat(spi_sampc); 
    spi_sort = sort(spi_samp);
    spi_res = sort(cell2mat(reshape(timernd,[],1)));

%% Any two spindle events less than 1 s, remove it!
     if length(spi_sort) > 0
        SPI_time = []; j=0;         
        for i = 1 : (length(spi_sort)-1)
            A = spi_sort(i+1); 
            B = spi_sort(i);
            tmp = spi_sort(i+1) - spi_sort(i);
            if tmp >= 0.2*Rs
                j = j+1;
                SPI_time(j) = spi_sort(i);
                SPI_RES(j) = spi_res(i);
            end 
        end 

        SPI_time = sort(SPI_time'); 
        SPI_RES = sort(SPI_RES'); 

     else
         SPI_time = []; SPI_RES = [];        
     end

%% If the trials are really huge, MATLAB goes crazy. Hence, divide trials into bunch of 800 trials.  
    
    if length(SPI_time) > 0 
        disp('SPI time > 0')        
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
        clear trial_start trial_end
 
        %Now for each of the small chunks of 800 trials, do PGD and store
        
        for qqq = 1 % length(trialS(:,1))      
        
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
            resample_factor=10;
            Data_hil_tmp_resample=resample (Data_hil_tmp', 1,resample_factor)';
            data_fil=reshape(Data_fil_tmp,CH_N,duration,trial_N);
            data_hil_all=reshape(Data_hil_tmp_resample,CH_N,duration/resample_factor,trial_N);

%             t=(1:duration/resample_factor)./(Rs/resample_factor)-0.5;
%             timezero = round(length(t+1)/2);
            midphiT = round(length(data_hil_all(1,:,1))/2);

    %% Check for bad channels INSIDE the cluster 
    
            indiv_clus = cluster_set(mm,:);  %does not need NGmap, already organized 
            clust_hil = data_hil_all(indiv_clus, :, :); %Now extracting the hilbert data for those cluster channels 
    
            kt=0;
            badindex = 0;
            for nn = 1:length(empsp)

              if ismember(empsp(nn),NGmap(indiv_clus)) == 1        %If bad channel falls inside this cluster, make a badindex
                  kt = kt+1;
                  disp(['Bad Channel inside clus=', num2str(mm)])
                  badindex(kt) = find(NGmap(indiv_clus) == empsp(nn));
              else          
                  nn = nn+1;
              end
            end 
 
            for trial = 1 : trial_N

                for ii = 1:length(indiv_clus)
                    instphase(ii,1) = unwrap((angle(clust_hil(ii,midphiT,trial)))); %Filling my good 3x3 cluster with phase values
                end

                if badindex > 0 
%                     
                    [cluster_set_interp, removerow, removecol]  = inpaintany(mm, gridsize);
                                       
                    indiv_clus_interp = cluster_set_interp;                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %New set of bad indices for 4x4 cluster                    
                    kt=0;
                    badindex_new = 0;
                    for nn = 1:length(empsp)

                      if ismember(empsp(nn),NGmap(indiv_clus_interp)) == 1        %If bad channel falls inside this cluster, make a badindex
                          kt = kt+1;
                          badindex_new(kt) = find(NGmap(indiv_clus_interp)== empsp(nn));
                      else          
                          nn = nn+1;
                      end
                    end              
                    
                    clust_hil_interp = data_hil_all(indiv_clus_interp, :, :); %taking Hilbert values in this new 4x4 cluster
                                       
                    for kk = 1:length(indiv_clus_interp)
                        instphase_interp(kk,1) = unwrap((angle(clust_hil_interp(kk,midphiT,trial)))); %fill this cluster with phase
                    end
                           
                    colg_i = gridsize+1; rowg_i = gridsize+1; 
                    instphase_interp([badindex_new],1) = NaN; %now Nan wherever there is a bad channel

                    gradgrid_old = reshape(instphase_interp,[colg_i,rowg_i])';

                    gradgrid = inpaint_nans(gradgrid_old,4); %inpaint this with method 4 
                    
                    gradgrid_final = gradgrid; 
                    gradgrid_final(removerow,:) = [];  %Now go back to 3x3 cluster
                    gradgrid_final(:,removecol) = []; 
                    
                    
                    INPAINT(trial).four = gradgrid_old; 
                    INPAINT(trial).paint = gradgrid; 
                    INPAINT(trial).final = gradgrid_final; 
                    
                    %FINAL PGD GRID BEFORE COMPUTATION
                    gradgrid_PGD = gradgrid_final; 
                     
                elseif badindex == 0
                    
                    %WITHOUT NANs, FINAL PGD GRID
                    gradgrid_PGD = reshape(instphase,[colg,rowg])';

                end
                
                [px, py] = gradient(gradgrid_PGD);

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

            end             

            Cluster_res(mm).clus = num2str(cluster_set(mm,:));
            Cluster_res(mm).RES((trial_start:trial_end),1) = SPI_RES(trial_start:trial_end);

            Cluster_res(mm).SPI((trial_start:trial_end),1) = SPI_time(trial_start:trial_end);
            Cluster_res(mm).PGD((trial_start:trial_end),1) = PGD;
            Cluster_res(mm).angles((trial_start:trial_end),1) = dirc;
            
            tmp1 = Cluster_res(mm).PGD;
            largePGD1 = find(tmp1>PGD_thresh); 

            percent = 100*(length(largePGD1)/length(tmp1));

            Cluster_res(mm).sigper = percent;
                
            clear PGD
            clear dirc
           
        end 
            
    elseif length(SPI_time) <=0 
        
        disp('SPI_time is less 0')                    
        Cluster_res(mm).clus = num2str(cluster_set(mm,:));
        Cluster_res(mm).RES = SPI_RES;
        Cluster_res(mm).SPI = SPI_time;
        Cluster_res(mm).PGD = [];
        Cluster_res(mm).angles = [];           
        Cluster_res(mm).sigper = [];
                
    end 
                
 
end
disp('PGD Calculation Done') 
%%

%% Plot one cluster
% close all 
% 
% % maxxr = 290; step = maxxr/2; 
% % 
%  dirc = Cluster_res(mm).angles;
% % direcsig = dirc(largePGD1);           %angle of direction vector for large PGDs 
% % 
% figure(1)
% polarhistogram(dirc,30) 
% pax = gca;
% % rticks([0:step:maxxr])
% thetaticks([0 90 180 270])
% pax.ThetaDir = 'clockwise';
% pax.ThetaZeroLocation = 'right';        
% % rlim([0  maxxr])
% % title('Cluster 3 4 5 (5 12 19 missing)')
% % 
% % print -depsc -tiff -r300 -painters Cluster3_Inpaint_5_12_19_missing.eps
% % 
% % 
% % 

%%


 save(savename,'Cluster_res')
    

     