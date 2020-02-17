%% Prawesh Dahal
% May 8, 2019 
%Detect spindles for lfp- wavelet - 

close all
clear all
clc 

%%

lfp_file = dir('*.lfp');
CH_N=xml2CH_N(cat(2,lfp_file.name(1:end-4),'.xml'));
[~, fbasename, ~] = fileparts(lfp_file.name);

spi_file = strcat(fbasename, '_SPI_RES'); 

%%
tic

parfor i = 1 : 128

[spindles] = Spindle_Detect_wavelet(CH_N, i)
 
spindle(i).res = spindles; 

spindle_events = strcat(fbasename,'_', num2str(i), '.spi.evt');
channelID = i-1; 
SaveRippleEvents(spindle_events, spindles, channelID);

mkdir SPI_Events
movefile([strcat(fbasename,'_', num2str(i), '.spi.evt')], strcat(pwd,'/','SPI_Events')) 

end 

toc
save(spi_file, 'spindle');


 