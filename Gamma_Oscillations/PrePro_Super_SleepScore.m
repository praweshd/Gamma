%% PrePro_Super_SleepScore
%Adapted from buzcode: SleepScoreMaster
%July 15, 2019 
%Prawesh Dahal (TNL)
%
%EMGCH option1
%[sessionInfo.AnatGrps(1).Channels(1:2) SWChannels(1:2)
%ThetaChannels(1:2)];

function PrePro_Super_SleepScore(ACCM, thetagroup, swgroup, emgroup, truelength)
 
warning('off','all')

basePath = pwd;
slash_index=find(basePath=='/');
currentFolder=basePath(slash_index(end)+1:end);

disp(['Auto Sleep Scoring for file ', currentFolder])

sessionInfo = LoadParameters(basePath);

ThetaChannels = sessionInfo.AnatGrps(thetagroup).Channels(1:end);
disp(['Picking HP CH: ', num2str(ThetaChannels) ])

SWChannels = sessionInfo.AnatGrps(swgroup).Channels(1:end);
disp(['Picking CX CH: ', num2str(SWChannels) ])

EMGCH = sessionInfo.AnatGrps(emgroup).Channels(1:3);
disp(['Picking EMG  : ', num2str(EMGCH) ])

SleepScoreMaster(basePath,EMGCH,ACCM,'SWChannels',SWChannels,'ThetaChannels',ThetaChannels)
filname = strcat(currentFolder, '.SleepState.states.mat');
 
load(filname);
states = SleepState.idx.states; 
states = states';  

if length(states) > truelength
    states = states(1,1:truelength);
elseif length(states) < truelength
    states(1, length(states):truelength) = 0; 

disp('Saving autostates .mat')

filename = dir('*.lfp');
fbasename = regexprep(filename.name,'.lfp','');

savename = [fbasename, '-autostates.mat']; 
save(savename,'states')

end

