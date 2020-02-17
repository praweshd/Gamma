%% Sleep_Score_Accelerometer
%Adapted from buzcode: SleepScoreMaster
%July 15, 2019 
%Prawesh Dahal (TNL)
%%
close all
clear all
clc
%%
% {Parameters}
%If this is manual sleep scored
state_mat = dir('*-teststates*');
load (state_mat.name);
StateIntervals = ConvertStatesVectorToIntervalSets(states);                 % 6 Intervalsets representing sleep states
REM = StateIntervals{5}; 
NREM = StateIntervals{3};
WAKE = StateIntervals{1};
state = NREM;

%Load LFP Filename
filename = dir('*.lfp');
fbasename = regexprep(filename.name,'.lfp','');

F_in = filename.name;
CH_N=xml2CH_N(cat(2, F_in(1:end-4),'.xml'));

truelength = length(states); 

%% Read accelerometer channel and resample
ACC_ch = 202; 
basePath = pwd;
slash_index=find(basePath=='/');
currentFolder=basePath(slash_index(end)+1:end);
sel_data = readmulti(strcat(currentFolder, '.lfp'), CH_N, ACC_ch); 
dLfp = sel_data(1:1250:end); 
%%
dLfp = dLfp(1:truelength);

%% High Pass Filter to get rid of drift
Rs = 1250;
n=3; Wn=0.1;
[b,a]=butter(n,2*Wn/Rs,'high');
fil_acc=filtfilt(b,a,dLfp);

%Rectify 
yy = ((fil_acc).^2);

%Take derivative
ydif = diff(yy);

%Rectify the derivative
drec = (ydif).^2;

%% Choose threshold from the test subplots
% Detect ripple periods by thresholding normalized squared signal
lowThresholdFactor = 1.8e+12; 
thresholded = drec > lowThresholdFactor;
start = find(diff(thresholded)>0);
stop = find(diff(thresholded)<0);

%% FIRST PASS
% Exclude last ripple if it is incomplete
if isempty(start),
	disp('Detection by thresholding failed');
    ripples = [0 0 0 0];
	return
else
    if length(stop) == length(start)-1,
        start = start(1:end-1);
    end
    % Exclude first ripple if it is incomplete
    if length(stop)-1 == length(start),
        stop = stop(2:end);
    end
    % Correct special case when both first and last ripples are incomplete
    if start(1) > stop(1),
        stop(1) = [];
        start(end) = [];
    end
    firstPass = [start,stop];
%     disp(['After detection by thresholding: ' num2str(length(firstPass)) ' events.']);
end
%% SECOND PASS

% Merge ripples if inter-ripple period is too short
%minInterRippleSamples = minInterRippleInterval/1000*frequency;
minInterRippleSamples = 40; 
secondPass = [];
ripple = firstPass(1,:);
for i = 2:size(firstPass,1)
	if firstPass(i,1) - ripple(2) < minInterRippleSamples,
		% Merge
		ripple = [ripple(1) firstPass(i,2)];
	else
		secondPass = [secondPass ; ripple];
		ripple = firstPass(i,:);
	end
end
secondPass = [secondPass ; ripple];
if isempty(secondPass),
	disp('Ripple merge failed');
    ripples = [0 0 0 0];
	return
else
% 	disp(['After ripple merge: ' num2str(length(secondPass)) ' events.']);
end
%% THIRD PASS
% Discard ripples with a peak power < highThresholdFactor
thirdPass = [];
peakNormalizedPower = [];
highThresholdFactor = 0.6*10^(13);
% disp(['peak threshold value: ' num2str(highThresholdFactor)]);

for i = 1:size(secondPass,1)
	[maxValue,maxIndex] = max(drec([secondPass(i,1):secondPass(i,2)]));
	if maxValue > highThresholdFactor
		thirdPass = [thirdPass ; secondPass(i,:)];
		peakNormalizedPower = [peakNormalizedPower ; maxValue];
	end
end

%%
wakestates = zeros(length(dLfp),1); 
for i = 1:length(thirdPass)    
    wakestates(thirdPass(i,1):thirdPass(i,2)) = 1;     
end 

%% Check relation of states and acc data 
close all
Rs = 1250; 
wini = 3000; 
winf = 4000; 
spnum = 6;

fig=figure_ctrl('Time trace',1500,1000);
SP1 = subaxis (spnum,1,1, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
plot(dLfp, 'b')
xlim([wini winf])

% SP3 = subaxis(spnum,1,3, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
% plot((fil_acc).^2, 'b'); 
% xlim([wini winf]) 
% 
% SP4 = subaxis(spnum,1,4, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
% plot(ydif,'b'); 
% xlim([wini winf]) 
  
SP2 = subaxis(spnum,1,2, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
plot(drec,'b'); 
xlim([wini winf]); ylim([0 6*10^(12)]); title('Fil->Rec->Diff->Rec');

SP3 = subaxis(spnum,1,3, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
plot(drec,'b'); hold on
plot(start, ones(length(start),1)*lowThresholdFactor,'g*')
plot(stop, ones(length(start),1)*lowThresholdFactor,'r*')
plot( ones(length(drec),1)*lowThresholdFactor,'k')
plot( ones(length(drec),1)*highThresholdFactor,'c')
xlim([wini winf]); ylim([0 6*10^(12)])

SP4 = subaxis(spnum,1,4, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
plot(drec,'b'); hold on
plot(secondPass(:,1), ones(length(secondPass(:,1)))*lowThresholdFactor,'g*')
plot(secondPass(:,2), ones(length(secondPass(:,2)))*lowThresholdFactor,'r*')
plot( ones(length(drec),1)*lowThresholdFactor,'k')
plot( ones(length(drec),1)*highThresholdFactor,'c')
xlim([wini winf]); ylim([0 2*10^(13)])

SP5= subaxis(spnum,1,5, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
plot(wakestates, 'k'); 
xlim([wini winf]); ylim([0 2]) 

SP6 = subaxis(spnum,1,6, 'Spacing', 0.02, 'Padding', 0.02, 'Margin', 0.02);
plot(states, 'b'); 
xlim([wini winf]); 

linkaxes([SP1 SP2 SP3 SP4 SP5 SP6], 'x')

%% Finally RUN sleepscoring
% 
ACCM = wakestates;
thetagroup = 11;
swgroup = 10;
emgroup = 10; 
  
PrePro_Super_SleepScore(ACCM, thetagroup, swgroup, emgroup, truelength)

%%
