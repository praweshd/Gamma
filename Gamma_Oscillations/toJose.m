%% Gamma Detection on Neurogrid 
% Prawesh Dahal
% June 17, 2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=============================================================
%=============================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Main Parameters Run This Always for all sections
close all
clear all
clc

%Load data from desired STATES
state_mat = dir('*-states*');
load (state_mat.name);
StateIntervals = ConvertStatesVectorToIntervalSets(states);                 % 6 Intervalsets representing sleep states
REM = StateIntervals{5}; 
NREM = StateIntervals{3};
WAKE = StateIntervals{1};
state = NREM;

lfp_file = dir('*.lfp');
CH_N=xml2CH_N(cat(2,lfp_file.name(1:end-4),'.xml'));

%Load Bad channels
badch = dir('*bad_CH*');
load(badch.name);
CH_noise = bad_ch'; 

%Load grid info 
NG = dir('NG2019_info.mat'); 
load(NG.name); 

NGmap = NG2019_info.map; 
NGdim = NG2019_info.dim;

NREM_dur = length(find(states==3)); 

num_CH = 205; %total num of channels
gam_CH = 10;  %desired sample channel 

Rs = 1250;

%% I - This section of code is if you already have detected start-end times of the event and want to 
%look around the window of that event
% Plot raw lfp, filtered and rectified sample traces
 
%Load detected event res file
detec_file = dir('*_all_gamma.mat');
load(detec_file.name); 


Rs = 1250;
gam = gammas(1).res; %gammas is my res struct file 
%time = round(gam(:,1)*Rs);

avg_res = (gam(:,1) + gam(:,3))/2;
time = round(avg_res*Rs);

duration = 0.5*Rs; 
t=(1:duration)./Rs-2.5;

durplot = 0.1*Rs; 
start_point=round((duration/2) - durplot );
end_point= round(start_point   + 2*durplot);

trial_N=400;
%Extract data
data = Dat_tracker(lfp_file.name,time(1:trial_N),duration,CH_N);

%Bad channels NaN them
data(CH_noise,:,:) = NaN;

Wn = [20 60]; resample_factor = 1; 

%Filter and Rectify the raw data
[data_fil, data_hil] = data_filhil(data,Rs,Wn,trial_N,duration,CH_N,resample_factor);

%Trigger Average
data_hil_avg=mean(abs(data_hil),3);
gam_hil=data_hil_avg(:,start_point:end_point);
gam_pwr=sum(gam_hil,2)./(length (start_point:end_point));   

%% Plotting all the three sets - RAW, FILTERED AND RECTIFIED

%Choose for what set of channels
CH_sel_all = [25:50];
trial = 1;

traceplot(data,trial,CH_sel_all)
title('Raw')

traceplot(data_fil,trial,CH_sel_all)
title('Filtered')

traceplot(abs(data_hil),trial,CH_sel_all)
title('Hilbert')

traceplot(data_hil_avg,trial,CH_sel_all)
title('Avg Hilbert all trials')

traceplot(gam_hil,trial,CH_sel_all)
title('Start-End')


%% II - This section is when you have not detected events yet and just want to
%familiarize with the raw lfp. See the effects of whitening, smoothening
%and look at the wavelet. 
close all
clc

%Load raw lfp
tic
[~, fbasename, ~] = fileparts(lfp_file.name);
lfp = LoadLfp(fbasename,CH_N,gam_CH); 
raw_lfp = [Range(Restrict(lfp,state),'s') Data(Restrict(lfp,state))];
toc

Nsamp = length(raw_lfp); 
resample_factor = 10; 
R_s =Rs/resample_factor;
%%
%Parameters
freq = 1:120; 

%Plotting window
wini = R_s; 
winf = 4*R_s; %the final plot you make is between 1s - 5s

%Smoothing window
smooth_win= round(R_s*(0.01));

arburg_n=2;

lowThresholdFactor = -50; % Ripple envoloppe must exceed lowThresholdFactor*stdev
highThresholdFactor = 0; %

gamb = [30 60];
lowb = [1 16] ;
highb = [80 120]; 

maxg = 300;
ming = 50;
interg = 100; 

%%
% Downsample
dLfpr = raw_lfp(:,2);
dLfp = resample(dLfpr, R_s, Rs); 
rg = raw_lfp(:,1);
rg = rg(1:resample_factor:end);
%%

%Whitened and Wavelet 
b = arburg(dLfp,arburg_n);
y = Filter0In(b, dLfp);
[SS,freq,~] = awt_freqlist(y,R_s,freq,'Gabor');
dd = abs(SS);

%%
%Smoothen  
S_smooth = smoothdata(dd,2,'gaussian',smooth_win);

%%
%Low Ratio 
temp = S_smooth;
gam_band = sum(temp(:,30:60), 2);
low_band = sum(temp(:, 1:16), 2);
high_band = sum(temp(:, 70:120), 2);
  
gam_ratio  =  (gam_band - low_band - high_band) ./ (gam_band + low_band + high_band);  

%%
%Plot: Time Trace
close all
clc
fig=figure_ctrl('Time trace',1500,1000);
SP1 = subaxis (4,1,1, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
plot(dLfp, 'b')
xlim([wini winf])
 

%  Plot: Whitened
SP2 = subaxis(4,1,2, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
imagesc(dd'); axis xy; colormap jet; caxis([0 80]);
xlim([wini winf]); ylim([0 80])

% Plot Smooth
SP3 = subaxis(4,1,3, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
imagesc(S_smooth'); axis xy; colormap jet; caxis([0 80]);
xlim([wini winf]); ylim([0 80])

%Plot Gam Ratio
SP4 = subaxis(4,1,4, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
plot(gam_ratio, 'b')
xlim([wini winf])
hold on
plot((1:Nsamp), ones(length((1:Nsamp)),1)*highThresholdFactor,'r')
%plot((1:Nsamp), ones(length((1:Nsamp)),1)*lowThresholdFactor,'g')

linkaxes([SP1 SP2 SP3 SP4], 'x')
 