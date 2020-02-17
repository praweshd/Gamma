
function traceplot(data, trial, CH_sel_all,  color)
%  
% %% %PGDs
% allpgd = dir('*PGD_NAN_M4*');
% load(allpgd.name);
% 
% filename = dir('*.lfp');                %Get filename  
% F_in = filename.name;
% CH_N=xml2CH_N(cat(2, F_in(1:end-4),'.xml'));
% 
% timespi =  25272.946 + 1.6 ; 
% %%
% Rs = 1250;
% num_CH = 139;           %input the number of channels
% CH_sel_all = [1:64];  %select the channels you want to display
% time = timespi;      % mid point of the interval that you want to display in seconds!
% frame_time =5;       % length of interval that you want to display in seconds!5
figure_name = 'barreltrace';  % how you want to label the figure!
neg_yaxis_val = 2500;   %adjust these based on the amplitude of your traces empirically
pos_yaxis_val = 2500;
% 
% filename = dir('*.lfp');
% [pathstr, fbasename, fileSuffix] = fileparts(filename.name);
% data = Dat_tracker(filename.name, round(time*Rs), frame_time*Rs,num_CH);
% 

if length(size(data)) == 3 
    
    data=data(CH_sel_all,:, trial);
    
elseif length(size(data)) == 2
    
    data=data(CH_sel_all,:);
    
end 
 


for ii = 1:length(CH_sel_all)
    plot (data(ii,:) -ii*300,color); 
    grid on
    axis_cleaner; axis tight
    hold on
end
% 
% % print('RAW', '-depsc', '-tiff','-r300', '-painters') 
% 
% %% Filter
% 
% % n=3;
% % Wn=[10 15];
% % [b,a]=butter(n,2*Wn/Rs,'bandpass'); 
% % data=filtfilt (b,a,data')';
% figure1 = figure_ctrl('Filtered spindle', 500, 1000);
% 
% 
% for ii = 1:length(CH_sel_all)
%     plot (data_fil(ii,:,trial) -ii*250,'k'); 
%     axis_cleaner; axis tight
%     hold on
% end
% 
% % print('Filtered', '-depsc', '-tiff','-r300', '-painters') 

% 
% %% Only plot
% 
% Rs = 1250;
% num_CH = 139;           %input the number of channels
% %CH_sel_few = [4:10 11 17 18 22 23 25 30 31 34:40 44 45 ];  %select the channels you want to display
% CH_sel_few = [25:32 ];
% time = time spi;        % mid point of the interval that you want to display in seconds!
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


end

