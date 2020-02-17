%% PLOT RAW TRACES

close all
clear all
clc

 %%
%Load Bad channels
badch = dir('*bad_CH*');
load(badch.name);
CH_noise = bad_ch'; 
%%
  
filename = dir('*.lfp');                %Get filename  
F_in = filename.name;
CH_N=xml2CH_N(cat(2, F_in(1:end-4),'.xml'));

timespi =  13209.138;
%%
close all
Rs = 1250;
CH_sel_all = [1:120];  %select the channels you want to display
time = timespi;      % mid point of the interval that you want to display in seconds!
frame_time =0.3;       % length of interval that you want to display in seconds!5
figure_name = 'barreltrace';  % how you want to label the figure!
neg_yaxis_val = 2500;   %adjust these based on the amplitude of your traces empirically
pos_yaxis_val = 2500;

filename = dir('*.lfp');
[pathstr, fbasename, fileSuffix] = fileparts(filename.name);
data = Dat_tracker(filename.name, round(time*Rs), frame_time*Rs,CH_N);

%%
data=data(CH_sel_all,:);

% n=3;
% Wn=[60 80];
% [b,a]=butter(n,2*Wn/Rs,'bandpass'); 
% data=filtfilt (b,a,data')';

figure1 = figure_ctrl('RAW', 300, 800);

k=0;
for ii = 1:length(CH_sel_all)
    if ~ismember(CH_sel_all(ii),bad_ch)
        k=k+1;
        plot (data(ii,:) -k*250,'k'); 
        axis_cleaner; axis tight
        hold on
    end
end

print(['traceR_', num2str(round(timespi))], '-depsc', '-tiff','-r300', '-painters') 


n=3;
Wn=[60 80];
[b,a]=butter(n,2*Wn/Rs,'bandpass'); 
datafil =filtfilt (b,a,data')';

figure2 = figure_ctrl('FIL', 300, 800);
k=0;
for ii = 1:length(CH_sel_all)
    if ~ismember(CH_sel_all(ii),bad_ch) 
        k = k+1;
        plot (datafil(ii,:) -k*250,'k'); 
        axis_cleaner; axis tight
        hold on
    end
end

print(['traceF_', num2str(round(timespi))], '-depsc', '-tiff','-r300', '-painters') 

% %% Filter
% 
% n=3;
% Wn=[10 15];
% [b,a]=butter(n,2*Wn/Rs,'bandpass'); 
% data=filtfilt (b,a,data')';
% figure1 = figure_ctrl('Filtered spindle', 500, 1000);
% 
% 
% for ii = 1:length(CH_sel_all)
%     plot (data(ii,:) -ii*250,'k'); 
%     axis_cleaner; axis tight
%     hold on
% end
% 
% print('Filtered', '-depsc', '-tiff','-r300', '-painters') 
% 
% 
% %% Only plot
% 
% Rs = 1250;
% num_CH = 139;           %input the number of channels
% %CH_sel_few = [4:10 11 17 18 22 23 25 30 31 34:40 44 45 ];  %select the channels you want to display
% CH_sel_few = [25:32 ];
% time = timespi;        % mid point of the interval that you want to display in seconds!
% frame_time = 1;       % length of interval that you want to display in seconds!5
% figure_name = 'barreltrace';  % how you want to label the figure!
% neg_yaxis_val = 2500;   %adjust these based on the amplitude of your traces empirically
% pos_yaxis_val = 2500;
% 
% filename = dir('*.lfp');
% [pathstr, fbasename, fileSuffix] = fileparts(filename.name);
% data = Dat_tracker(filename.name, round(time*Rs), frame_time*Rs,num_CH);
% data=data(CH_sel_few,:);
% 
% n=3;
% Wn=[10 15];
% [b,a]=butter(n,2*Wn/Rs,'bandpass'); 
% data=filtfilt (b,a,data')';
% 
% % 
% % n=3;
% % Wn=2;
% % [b,a]=butter(n,2*Wn/Rs,'high');
% % data=filtfilt(b,a,data')';
% figure3 = figure_ctrl('sample spindle', 500, 1000);
% 
% for ii = 1:length(CH_sel_few)
%     plot (data(ii,:) -ii*250,'k'); 
%     axis_cleaner; axis tight
%     hold on
% end
% 
% 
% title('1s trace')
% print('1s trace', '-depsc', '-tiff','-r300', '-painters') 
% 
% 
