%% Gamma Detection on Neurogrid 
% Prawesh Dahal
% June 17, 2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=============================================================
%=============================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Main Parameters Run This 
close all
clear all
clc

state_mat = dir('*-states*');
load (state_mat.name);
StateIntervals = ConvertStatesVectorToIntervalSets(states);                 % 6 Intervalsets representing sleep states
REM = StateIntervals{5}; 
NREM = StateIntervals{3};
WAKE = StateIntervals{1};
state = NREM;

%Parameters
Rs = 1250; 
resample_factor = 1; 
R_s =Rs/resample_factor;
freq = 1:350; 

wini = R_s; 
winf = 1*R_s; 
smooth_win= round(R_s*(0.01));

arburg_n=2;

% lowThresholdFactor = -50; % Ripple envoloppe must exceed lowThresholdFactor*stdev
% highThresholdFactor = 0; %
% 
% gamb = [60 100];
% lowb = [1 40] ;
% highb = [120 140]; 
% 
% maxg = 300;
% ming = 50;
% interg = 100; 

tic
lfp_file = dir('*.lfp');
CH_N=xml2CH_N(cat(2,lfp_file.name(1:end-4),'.xml'));


%% ================================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================================
% 1. WAVELET BASED - DOWNSAMPLING
%% Run this if downsampling
tic
[~, fbasename, ~] = fileparts(lfp_file.name);
lfp = LoadLfp(fbasename,CH_N,gam_CH); 
raw_lfp = [Range(Restrict(lfp,state),'s') Data(Restrict(lfp,state))];
toc
%%

Nsamp = length(raw_lfp); 
resample_factor = 10; 
R_s =Rs/resample_factor;

% Downsample
dLfpr = raw_lfp(:,2);
dLfp = resample(dLfpr, R_s, Rs); 
rg = raw_lfp(:,1);
rg = rg(1:resample_factor:end);

dLfp = dLfpr;

%Whitened and Wavelet 
b = arburg(dLfp,arburg_n);
y = Filter0In(b, dLfp);
[SS,freq,~] = awt_freqlist(y,R_s,freq,'Gabor');
dd = abs(SS);

%Smoothen  
S_smooth = smoothdata(dd,2,'gaussian',smooth_win);

%Low Ratio 
temp = S_smooth;
gam_band = sum(temp(:,30:60), 2);
low_band = sum(temp(:, 1:16), 2);
high_band = sum(temp(:, 70:120), 2);
  
gam_ratio  = gam_band ;%(gam_band - low_band - high_band) ./ (gam_band + low_band + high_band);  

%Plot: Time Trace
close all
clc
fig=figure_ctrl('Time trace',1500,1000);
SP1 = subaxis (5,1,1, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
plot(dLfp, 'b')
xlim([wini winf])

% Plot: Wavelet
SP2 = subaxis(5,1,2, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
imagesc(d'); axis xy; colormap jet; caxis([0 80]);
xlim([wini winf]); ylim([0 80])

%  Plot: Whitened
SP3 = subaxis(5,1,3, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
imagesc(dd'); axis xy; colormap jet; caxis([0 80]);
xlim([wini winf]); ylim([0 80])

% Plot Smooth
SP4 = subaxis(5,1,4, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
imagesc(S_smooth'); axis xy; colormap jet; caxis([0 80]);
xlim([wini winf]); ylim([0 80])

%Plot Gam Ratio
SP5 = subaxis(5,1,5, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
plot(gam_ratio, 'b')
xlim([wini winf])
hold on
plot((1:Nsamp), ones(length((1:Nsamp)),1)*highThresholdFactor,'r')
plot((1:Nsamp), ones(length((1:Nsamp)),1)*lowThresholdFactor,'g')

linkaxes([SP1 SP2 SP3 SP4 SP5], 'x')


%% ================================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================================
%% LFP TO LOW/HIGH GAMMA RATIO FILE 
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this is the .spi equivalent 
[GAM] = dat2gam(dLfpr,resample_factor,arburg_n,Rs,freq,smooth_win,gamb,lowb,highb);
lFpgam = [rg GAM'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamb = [60 100];
lowb = [1 16] ;
highb = [110 120];

[GAMH] = dat2gam(dLfpr,resample_factor,arburg_n,Rs,freq,smooth_win,gamb,lowb,highb);
lFpgam = [rg GAMH'];
maxg = 300;
ming = 20; 
interg = 100; 

% Event Detection
[gam_evt] = FindGamma_wav(lFpgam, 'thresholds', [lowThresholdFactor highThresholdFactor], 'durations', [interg maxg ming], 'frequency', 125); 

gamma_file = strcat(fbasename, '_', num2str(gam_CH), 'gamma_wav'); 
save(gamma_file, 'gam_evt');
gamma_events = strcat(fbasename,'_', num2str(gam_CH), '.gam.evt');
channelID = gam_CH-1; 
SaveRippleEvents(gamma_events, gam_evt, channelID);

disp('Low Gamma Event file has been saved')


%% ================================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%================================================================================================== 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================================
% 1. WAVELET BASED - NO DOWNSAMPLING 
%% Spi Lfp compare 
% lfp2gam_detect_wavelet(freq,gamb,lowb,highb); 
 
close all
clc 

lowThresholdFactor = -50; % Ripple envoloppe must exceed lowThresholdFactor*stdev
highThresholdFactor = 0; %

num_CH = 1; gam_CH = 1; 
spi_file = dir('*.spi');
[~, fbasename, ~] = fileparts(spi_file.name);

 
% 1725.038 4825.797 4890.524 4899.957
time = round(1729.576*Rs);

duration=5*Rs;

dataLFP=Dat_tracker(lfp_file.name,time,duration,CH_N);
dataLFP=squeeze(dataLFP(gam_CH,:));

dataSPI = Dat_tracker(spi_file.name,time,duration,1);
dataSPI=squeeze(dataSPI(1,:));
 
wN = length(dataLFP); 

b = arburg(dataLFP,arburg_n);
y = Filter0In(b, dataLFP);

[SS,freq,~] = awt_freqlist(y,R_s,freq,'Gabor');
dd = abs(SS);

smooth_win= round(R_s*(0.01));
S_smooth = smoothdata(dd,2,'gaussian',smooth_win);

%time= round(2731.045*Rs);
%Plot: Time Trace? 
close all
clc
fig=figure_ctrl('Time trace no downsample',1500,1000);
SP1 = subaxis (3,1,1, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
plot(dataLFP, 'b')
axis tight

SP2 = subaxis (3,1,2, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
imagesc(dd'); axis xy; colormap jet; %caxis([0 15]); 

SP3 = subaxis (3,1,3, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
plot(dataSPI, 'b')
hold on
plot((1:wN), ones(length((1:wN)),1)*highThresholdFactor,'r')
plot((1:wN), ones(length((1:wN)),1)*lowThresholdFactor,'g')
axis tight
 
%% If above looks right - Gamma Event Detect
 
maxg = 300;
ming = 50;
interg = 100; 

gam_CH = 10; 

% Load file and restrict to state
filename = dir('*.spi');
[~, fbasename, ~] = fileparts(filename.name);
lfp = LoadGamma(fbasename,'.spi',1,1, 'frequency', 1250); 
restricted_lfp = [Range(Restrict(lfp,state),'s') Data(Restrict(lfp,state))];

% Event Detection
[gam_evt] = FindGamma_wav(restricted_lfp, 'thresholds', [lowThresholdFactor highThresholdFactor], 'durations', [interg maxg ming], 'frequency', 1250); 

gamma_file = strcat(fbasename, '_W', num2str(gam_CH), 'gamma_wav'); 
save(gamma_file, 'gam_evt');
gamma_events = strcat(fbasename,'_W', num2str(gam_CH), '.gam.evt');
channelID = gam_CH-1; 
SaveRippleEvents(gamma_events, gam_evt, channelID);

disp('Low Gamma Event file has been saved')

%% ================================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================================
%% Run this for WITHOUT WAVELET - NO DOWNSAMPLING

close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%% GET STD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gam_CH = 18; %matlab channel

% Load entire file to first get the std
tic
[~, fbasename, ~] = fileparts(lfp_file.name);
lfp = LoadLfp(fbasename,CH_N,gam_CH); 
raw_lfp = [Range(Restrict(lfp,state),'s') Data(Restrict(lfp,state))];
toc

signal_lfp = raw_lfp(:,2);

Wn = [60 80]; % gamma freq range [Hz]
[b,a] = butter(3, 2*Wn/Rs, 'bandpass');
fil_sleep = filtfilt(b,a,signal_lfp);
sd = []; keep = [];
squaredSignal = fil_sleep.^2;
std_sig = squaredSignal(1:5*60*Rs);
[~,sd] = unity(std_sig,sd,keep);
%%
duration=1*Rs;

lowThresholdFactor = 1; % Ripple envoloppe must exceed lowThresholdFactor*stdev
highThresholdFactor = 3; %

% 1725.038 4825.797 4890.524 4899.957 1729.576 1732.270

%high gamma
% 1734.556 1868.241 1869.178 1871.286
% time = round(3929.490*Rs); %MAIN

time = round(482.446*Rs); 

dataLFP=Dat_tracker(lfp_file.name,time,duration,CH_N);
dataLFP=squeeze(dataLFP(gam_CH,:));
wN = length(dataLFP); 

b = arburg(dataLFP,arburg_n);
y = Filter0In(b, dataLFP);

[SS,freq,~] = awt_freqlist(y,R_s,freq,'Gabor');
dd = abs(SS);

signal_lfp = dataLFP;

% Bandpass Filter
Wn = [60 80]; % spindle freq range [Hz]
[b,a] = butter(3, 2*Wn/Rs, 'bandpass');
fil_sleep = filtfilt(b,a,signal_lfp);

% Square and normalize signal
 
squaredSignal = fil_sleep.^2;
[normalizedSquaredSignal,sd] = unity(squaredSignal,sd,keep);
 
close all
% %% Final Flow Plot 
fig=figure_ctrl('Detect Step',1500,1000);
picnum = 4; 
subaxis (picnum,1,1, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
plot( dataLFP ) ; axis tight
title('Without wavelet ratio')

subaxis (picnum,1,2, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
imagesc(dd'); axis xy; colormap jet; ylim([0 120]); caxis([0 8 ]);  

subaxis (picnum,1,3, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
plot( fil_sleep ); axis tight
 
subaxis (picnum,1,4, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
plot( normalizedSquaredSignal ); axis tight
hold on 
plot( ones(length(squaredSignal),1)*lowThresholdFactor,'g')
plot( ones(length(squaredSignal),1)*highThresholdFactor,'r')


 
% print([fbasename, '_gamspi'], '-depsc', '-tiff','-r300', '-painters')

%%  If above looks right - Gamma Event Detect
clc
maxg = 300;
ming = 20;
interg = 50; 

 
% Event Detection
[gam_evt] = FindGamma_I(raw_lfp,'short','thresholds', [lowThresholdFactor highThresholdFactor], 'durations', [interg maxg ming], 'frequency', 1250); 

gamma_file = strcat(fbasename, '_highGAM_', num2str(gam_CH), 'gamma_wav'); 
save(gamma_file, 'gam_evt');
gamma_events = strcat(fbasename,'_highGAM_', num2str(gam_CH), '.gam.evt');
channelID = gam_CH-1; 
SaveRippleEvents(gamma_events, gam_evt, channelID);
% 
% disp('No Wavelet Low Gamma Event file has been saved')



%% ===========================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================================

%% GRAND GAMMA DETECTION 

close all
clear all
clc

lfp_file = dir('*.lfp');
CH_N=xml2CH_N(cat(2,lfp_file.name(1:end-4),'.xml'));

tic
for i = 1 : CH_N 

[gamma] = Gamma_Detect(205, i); 

gammas(i).GAM_time = gamma; 

end 

toc

%%

for i = 1 : 5
    
    tmp(i) = gammas(i).GAM_time(1,1);
    
end 

%%

plot(tmp)
ylim([1700 1720])


%% CROSS-CORRELOGRAM 

ts1 = gammas(1).GAM_time(1:200);
ts2 = gammas(2).GAM_time(1:200);
window = [-0.3 0.3]; 

[tsOffsets, ts1idx, ts2idx] = crosscorrelogram(ts1, ts2, window);

hist(tsOffsets, 50)
% %% Downsample 
% 
% 
% 
% 
% %% Wavelet (no smoothening or whitening) 
% %  
% %     x=squeeze(data(CH,:,i));
% %    [S(:,:,i),freq,~] = awt_freqlist(x,Rs,freq,'Gabor');
% %     
% % end
% % fig=figure_ctrl(' SPI___S',800,1000);
% % for i=1:trial_N
% %     subaxis (8,1,i, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
% %      d = abs(squeeze(S(:,:,i)));
% %     imagesc(d'); axis xy; colormap jet; caxis([0 80]);
% % end
% % 

% 
% %% %% Downsample 
% close all
% clc
% Fs = 1/125; 
% 
% % Define thresholds
% lowThresholdFactor = 0.2;
% highThresholdFactor = 0.5;
% 
% 
% %     'durations'   min inter-ripple interval, max ripple duration and min ripple duration, in ms
% %                   (default = [30 100 20])
% 
% minInterRippleInterval = 30; % in ms
% maxRippleDuration = 300; % in ms
% minRippleDuration = 50; % in ms
% 
% dLfpr = restricted_lfp(:,2);
% dLfp = resample(dLfpr, 125, 1250); 
% rg = restricted_lfp(:,1);
% rg = rg(1:10:end);
% 
% time = 0:Fs:(length(dLfp)-1)*Fs;
% %time = rg;
%  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig = figure_ctrl('Raw LFP',1500,1000);
% plot(time, dLfp); title('Step 1: Downsampled to 125 Hz'); xlabel('time(s)')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % nfft = 1024;
% % fs = 125; 
% % window = nfft; 
% % noverlap = nfft/2;
% % 
% % [s,f,t] = spectrogram(dLfp,window,noverlap,[],fs);
% % 
% % % figure
% % % spectrogram(dLfp,window,noverlap,[],fs,'yaxis')
% % % colormap jet
% 
% %% Whiten
% % fig=figure_ctrl('Whitened',1500,1000);
% A = arburg(dLfp,2);
% Wdlfp = Filter0(A, dLfp);
% % plot(time, Wdlfp) 
% 
% %% Spectrogram  
% % [s,f,t] = spectrogram(Wdlfp,window,noverlap,[],fs);
% % 
% % figure
% % spectrogram(Wdlfp,window,noverlap,[],fs,'yaxis')
% % colormap jet

% 
% % Plot BP signal
% Rs = 125;
% time = 0:Fs:(length(fil_sleep(:,2))-1)*Fs;
% FIL =  fil_sleep(:,2);
% zoomup_i = 0.5*10^6; 
% zoomup_f = zoomup_i + 10*Rs; 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fig=figure_ctrl('Bandpassed',1500,1000);
% % subplot(2,1,1); plot(time, FIL); title(['Step 1: BP ', num2str(Wn(1)),'-', num2str(Wn(2)), ' Hz NREM signal']); xlabel('time (s)');
% % subplot(2,1,2); plot(time(zoomup_i:zoomup_f), FIL(zoomup_i:zoomup_f)); title('Zoom'); xlabel('time (s)'); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % [s,f,t] = spectrogram(FIL,window,noverlap,[],fs);
% % 
% % figure
% % spectrogram(FIL,window,noverlap,[],fs,'yaxis')
% % colormap jet
% %% Square
%   
% FILSQ = FIL.^2;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fig=figure_ctrl('Squared',1500,1000);
% % plot(time(zoomup_i:zoomup_f), FILSQ(zoomup_i:zoomup_f)); title('Step 2: Squared'); xlabel('time (s)'); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% Window
% windowLength = round(125/1250*11);
% window = ones(windowLength,1)/windowLength;
% 
% FILWIN = Filter0(window,sum(FILSQ,2)); 
% % 
% % fig=figure_ctrl('Window',1500,1000);
% % plot(time(zoomup_i:zoomup_f), FILWIN(zoomup_i:zoomup_f))
% % hold on
% % plot(time(zoomup_i:zoomup_f), ones(length(time(zoomup_i:zoomup_f)),1)*mean(FILWIN),'r')
% % title('Step 3: Window'); xlabel('time (s)'); 
% 

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fig=figure_ctrl('Normalize',1500,1000);
% % subplot(2,1,1)
% % plot(time, normalizedSquaredSignal)
% % hold on
% % plot(time, ones(length(time),1)*mean(normalizedSquaredSignal), 'r') 
% % ylim([-0.25 15]);
% % 
% % subplot(2,1,2)
% % plot(time(zoomup_i:zoomup_f), normalizedSquaredSignal(zoomup_i:zoomup_f))
% % hold on 
% % plot(time(zoomup_i:zoomup_f), ones(length(time(zoomup_i:zoomup_f)),1)*mean(normalizedSquaredSignal),'r')
% % title('Step 4: Normalized signal - mean = 0, std = 1'); 
% % xlabel('time (s)'); xlim([4000 4010]); ylim([-0.25 2]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% %% Final Flow Plot 
% fig=figure_ctrl('Detect Step',1500,1000);
% picnum = 4; 
% subplot(picnum,1,1)
% plot(time(zoomup_i:zoomup_f), dLfp(zoomup_i:zoomup_f))
% 
% subplot(picnum,1,2)
% plot(time(zoomup_i:zoomup_f), FIL(zoomup_i:zoomup_f))
% 
% subplot(picnum,1,3)
% plot(time(zoomup_i:zoomup_f), FILSQ(zoomup_i:zoomup_f))
% 
% subplot(picnum,1,4)
% plot(time(zoomup_i:zoomup_f), normalizedSquaredSignal(zoomup_i:zoomup_f))
% hold on 
% plot(time(zoomup_i:zoomup_f), ones(length(time(zoomup_i:zoomup_f)),1)*mean(normalizedSquaredSignal),'r')
% plot(time(zoomup_i:zoomup_f), ones(length(time(zoomup_i:zoomup_f)),1)*lowThresholdFactor,'g')
% plot(time(zoomup_i:zoomup_f), ones(length(time(zoomup_i:zoomup_f)),1)*highThresholdFactor,'c')
% xlim([4000 4010]); ylim([-0.25 2]);
% 
% %% Threshold (FIRST-PASS)
%  
% thresholded = normalizedSquaredSignal > lowThresholdFactor;
% start = find(diff(thresholded)>0);
% stop = find(diff(thresholded)<0);
% 
% % Exclude last ripple if it is incomplete
% if length(stop) == length(start)-1,
% 	start = start(1:end-1);
% end
% % Exclude first ripple if it is incomplete
% if length(stop)-1 == length(start),
%     stop = stop(2:end);
% end
% % Correct special case when both first and last ripples are incomplete
% if start(1) > stop(1),
% 	stop(1) = [];
% 	start(end) = [];
% end
% 
% firstPass = [start, stop];
% %% FIRST PASS PLOT 
% 
% start1 = firstPass(:,1)/Rs;
% stop1 = firstPass(:,2)/Rs;
% 
% fig=figure_ctrl('SPI Firstpass',1500,1000);
% plot(time(zoomup_i:zoomup_f), normalizedSquaredSignal(zoomup_i:zoomup_f))
% hold on 
% plot(time(zoomup_i:zoomup_f), ones(length(time(zoomup_i:zoomup_f)),1)*lowThresholdFactor,'k')
% plot(time(zoomup_i:zoomup_f), ones(length(time(zoomup_i:zoomup_f)),1)*highThresholdFactor,'c')
% plot(start1, ones(length(start1),1)*lowThresholdFactor,'g*')
% plot(stop1, ones(length(stop1),1)*lowThresholdFactor,'r*')
% xlim([4000 4010]); ylim([-0.1 2]); title('Firstpass')
% 
% 
% %% SECOND PASS
% minInterRippleInterval = 30; 
% 
% % Merge ripples if inter-ripple period is too short
% minInterRippleSamples = minInterRippleInterval/1000*Rs;
% secondPass = [];
% ripple = firstPass(1,:);
% for i = 2:size(firstPass,1)
% 	if firstPass(i,1) - ripple(2) < minInterRippleSamples,
% 		% Merge
% 		ripple = [ripple(1) firstPass(i,2)];
% 	else
% 		secondPass = [secondPass ; ripple];
% 		ripple = firstPass(i,:);
% 	end
% end
% secondPass = [secondPass ; ripple];
% 
% %% SECOND PASS PLOT
% 
% start2 = secondPass(:,1)/Rs;
% stop2 = secondPass(:,2)/Rs;
% 
% fig=figure_ctrl('SPI SecondPass',1500,1000);
% plot(time(zoomup_i:zoomup_f), normalizedSquaredSignal(zoomup_i:zoomup_f))
% hold on 
% plot(time(zoomup_i:zoomup_f), ones(length(time(zoomup_i:zoomup_f)),1)*lowThresholdFactor,'k')
% plot(time(zoomup_i:zoomup_f), ones(length(time(zoomup_i:zoomup_f)),1)*highThresholdFactor,'c')
% plot(start2, ones(length(start2),1)*lowThresholdFactor,'g*')
% plot(stop2, ones(length(stop2),1)*lowThresholdFactor,'r*')
% xlim([4000 4010]); ylim([-0.1 2]); title(' SecondPass')
% 
% %% THIRD PASS
%  
% % Discard ripples with a peak power < highThresholdFactor
% thirdPass = [];
% peakNormalizedPower = [];
% for i = 1:size(secondPass,1)
% 	[maxValue,maxIndex] = max(normalizedSquaredSignal([secondPass(i,1):secondPass(i,2)]));
% 	if maxValue > highThresholdFactor,
% 		thirdPass = [thirdPass ; secondPass(i,:)];
% 		peakNormalizedPower = [peakNormalizedPower ; maxValue];
% 	end
% end
% 
% 
% %% THIRD PASS PLOT 
% 
% start3 = thirdPass(:,1)/Rs;
% stop3 = thirdPass(:,2)/Rs;
% 
% fig=figure_ctrl('SPI ThirdPass',1500,1000);
% plot(time(zoomup_i:zoomup_f), normalizedSquaredSignal(zoomup_i:zoomup_f))
% hold on 
% plot(time(zoomup_i:zoomup_f), ones(length(time(zoomup_i:zoomup_f)),1)*lowThresholdFactor,'k')
% plot(time(zoomup_i:zoomup_f), ones(length(time(zoomup_i:zoomup_f)),1)*highThresholdFactor,'c')
% plot(start3, ones(length(start3),1)*lowThresholdFactor,'g*')
% plot(stop3, ones(length(stop3),1)*lowThresholdFactor,'r*')
% xlim([4000 4010]); ylim([-0.1 2]);title('  ThirdPass')
% 
% %% 
% % 
% % % Detect negative peak position for each ripple
% % peakPosition = zeros(size(thirdPass,1),1);
% % for i=1:size(thirdPass,1),
% % 	[minValue,minIndex] = min(signal(thirdPass(i,1):thirdPass(i,2)));
% % 	peakPosition(i) = minIndex + thirdPass(i,1) - 1;
% % end
% % 
% % % Discard ripples that are way too long
% % time = filtered(:,1);
% % ripples = [time(thirdPass(:,1)) time(peakPosition) time(thirdPass(:,2)) peakNormalizedPower];
% % duration = ripples(:,3)-ripples(:,1);
% % ripples(duration>maxRippleDuration/1000,:) = [];
% % disp(['After max duration test: ' num2str(size(ripples,1)) ' events.']);
% % 
% % %Discard ripples that are way too short
% % duration = ripples(:,3)-ripples(:,1);
% % ripples(duration<minRippleDuration/1000,:) = [];
% % disp(['After min duration test: ' num2str(size(ripples,1)) ' events.']);

%%


function y = Filter0(b,x)

if size(x,1) == 1,
	x = x(:);
end

if mod(length(b),2) ~= 1,
	error('filter order should be odd');
end

shift = (length(b)-1)/2;

[y0 z] = filter(b,1,x);

y = [y0(shift+1:end,:) ; z(1:shift,:)];

end 

function [U,stdA] = unity(A,sd,restrict)

if ~isempty(restrict),
	meanA = mean(A(restrict));
	stdA = std(A(restrict));
else
	meanA = mean(A);
	stdA = std(A);
end
if ~isempty(sd),
	stdA = sd;
end

U = (A - meanA)/stdA;

end 



